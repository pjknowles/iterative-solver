#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/RSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace xspace {
namespace detail {

template <class RS, class QS>
auto generate_candidates(const RS& rs, const QS& qs) {
  auto candidates = std::map<size_t, std::vector<size_t>>{};
  for (auto root : rs.working_set()) {
    candidates[root] = qs.modification_candidates(root);
  }
  return candidates;
}

template <class RS, class QS>
auto generate_candidates_flat_list(const RS& rs, const QS& qs) {
  auto candidates = std::vector<size_t>{};
  for (auto root : rs.working_set()) {
    auto croot = qs.modification_candidates(root);
    std::copy(begin(croot), end(croot), std::back_inserter(candidates));
  }
  std::sort(begin(candidates), end(candidates));
  return candidates;
}

// FIXME for now I removed merger pathway.
auto generate_pairs(const std::map<size_t, std::vector<size_t>>& candidates) {
  auto pairs = std::vector<std::pair<size_t, size_t>>{};
  for (const auto& root_candidates : candidates) {
    for (const auto& c : root_candidates.second) {
      pairs.emplace_back(c, c);
    }
    if (false) {
      const auto& c = root_candidates.second;
      if (c.size() == 1)
        pairs.emplace_back(c.front(), c.front());
      else
        for (size_t i = 1; i < c.size(); ++i)
          pairs.emplace_back(c[i - 1], c[i]);
    }
  }
  return pairs;
}

bool empty_candidates(const std::map<size_t, std::vector<size_t>>& candidates) {
  bool empty = true;
  for (const auto& c : candidates) {
    empty &= c.second.empty();
  }
  return empty;
}

/*!
 * @brief Performs Gram-Schmidt orthogonalisation
 *
 * For a vector set {v} with overlap matrix \f$ S_{ij}= v_i . v_j \f$, generate a linear transformation L to a new
 * orthogonal vector set \f$ u_i = \sum_j L_{ij} v_j \f$.
 *
 * u_0 = v_0
 * u_1 = v_1 - <v_1, u_0> / <u_0, u_0>
 * u_i = v_i -\sum_j=0^{i-1} <v_i, u_j> / <u_j, u_j> u_j
 *     = v_i -\sum_j=0^{i-1} <v_i, u_j> / <u_j, u_j> \sum_{k=0}^j L_{jk} v_k
 *     = v_i -\sum_j=0^{i-1} (\sum_{k=0}^{j} W_{ij} / N_{jj} L_{jk}) v_k
 *
 * Let the overlap between u and v be, \f$ W_{ij} = <v_i, u_j> \f$
 *
 * W_{ij} = <v_i, u_j> = \sum_{k=0}^{j} S_{ik} L_{jk}
 * N_{ii} = <u_i, u_i>
 *        = \sum_{k} L_{ik} <v_k, u_i>
 *        = \sum_{km} L_{ik} L_{im} S_{km}
 *
 * L_{0,j} = \delta_{i,j}
 * L_{i,k} = \delta_{i,k} - (1 + \delta_{i,k})(\sum_{j=0}^{k} W_{ij} / N_{jj} L_{jk})
 *
 * @param s overlap matrix
 * @param l row matrix with linear transformation into orthonormal basis. Upper triangular component is zero.
 * @returns norms of transformed vectors
 */
template <typename T>
std::vector<T> gram_schmidt(const Matrix<T>& s, Matrix<T>& l) {
  assert(s.rows() == s.cols());
  auto n = s.rows();
  l.fill(0);
  l.resize({n, n});
  auto norm = std::vector<double>{};
  auto w = std::vector<double>{};
  for (size_t i = 0; i < n; ++i) {
    w.assign(i, 0.);
    for (size_t j = 0; j < i; ++j) {
      for (size_t k = 0; k <= j; ++k) {
        w[j] += s(i, k) * l(j, k);
      }
    }
    for (size_t j = 0; j < i; ++j) {
      for (size_t k = 0; k <= j; ++k) {
        l(i, k) -= w[j] / norm[j] * l(j, k);
      }
    }
    l(i, i) = 1.;
    norm.emplace_back(0);
    for (size_t j = 0; j <= i; ++j) {
      for (size_t k = 0; k < j; ++k) {
        norm.back() += 2 * l(i, j) * l(i, k) * s(j, k);
      }
      norm.back() += l(i, j) * l(i, j) * s(j, j);
    }
  }
  std::transform(begin(norm), end(norm), begin(norm), [](auto el) { return std::sqrt(std::abs(el)); });
  return norm;
}

} // namespace detail

template <class R, class P, class Q, class ST>
void check_conditioning(XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>, ST>& xs, RSpace<R, Q, P>& rs,
                        QSpace<R, Q, P>& qs, PSpace<R, P>& ps, double svd_threshold, double norm_threshold,
                        Logger& logger) {
  logger.msg("xspace::check_conditioning", Logger::Trace);
  bool stable = false;
  auto candidates = detail::generate_candidates(rs, qs);
  while (!stable && !detail::empty_candidates(candidates)) {
    auto& s = xs.data[EqnData::S];
    auto svd = svd_system(xs.size(), array::Span<double>{&s(0, 0), s.size()}, svd_threshold);
    stable = svd.empty();
    if (!svd.empty()) {
      logger.msg("singular value = " + Logger::scientific(svd.front().value), Logger::Debug);
      if (logger.data_dump) {
        logger.msg("singular vector = ", begin(svd.front().v), end(svd.front().v), Logger::Info);
      }
      auto pairs = detail::generate_pairs(candidates);
      auto norms = std::vector<double>{};
      for (const auto& p : pairs) {
        const auto oQ = xs.dimensions().oQ;
        const auto& v = svd.front().v;
        const auto norm = v.at(oQ + p.first) * v.at(oQ + p.first) + v.at(oQ + p.second) * v.at(oQ + p.second);
        norms.emplace_back(norm);
      }
      auto it_max = std::max_element(begin(norms), end(norms));
      auto i_max = std::distance(begin(norms), it_max);
      stable = *it_max < norm_threshold;
      if (*it_max > norm_threshold) {
        qs.merge(pairs.at(i_max));
        xs.build_subspace(rs, qs, ps);
        candidates = detail::generate_candidates(rs, qs);
      }
    }
  }
}

/*!
 * @brief Generates linear transformation of X space to an orthogonal subspace and removes any redundant Q vectors.
 * @param xs X space container
 * @param rs R space container
 * @param qs Q space container
 * @param ps P space container
 * @param lin_trans  Linear transformation matrix to an orthogonal subspace
 * @param norm_threshold
 * @param logger
 */
template <class R, class P, class Q, class ST>
void check_conditioning_gram_schmidt(XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>, ST>& xs,
                                     RSpace<R, Q, P>& rs, QSpace<R, Q, P>& qs, PSpace<R, P>& ps, Matrix<ST>& lin_trans,
                                     double norm_threshold, Logger& logger) {
  logger.msg("xspace::check_conditioning_gram_schmidt", Logger::Trace);
  bool stable = false;
  auto candidates = detail::generate_candidates_flat_list(rs, qs);
  while (!stable && !candidates.empty()) {
    auto& s = xs.data[EqnData::S];
    auto norm = detail::gram_schmidt(s, lin_trans);
    if (logger.data_dump)
      logger.msg("norm after Gram-Schmidt = ", begin(norm), end(norm), Logger::Info);
    auto it_min = std::find_if(begin(norm), end(norm), [norm_threshold](auto el) { return el < norm_threshold; });
    size_t imin = norm.size();
    for (auto c : candidates) {
      if (norm[c] < norm_threshold) {
        imin = c;
        break;
      }
    }
    stable = (imin == norm.size());
    if (!stable) {
      logger.msg("removing candidate from q space i =" + std::to_string(imin) +
                     ", norm = " + Logger::scientific(norm[imin]),
                 Logger::Debug);
      qs.erase(imin);
      xs.build_subspace(rs, qs, ps);
      candidates = detail::generate_candidates_flat_list(rs, qs);
    }
  }
}

} // namespace xspace
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
