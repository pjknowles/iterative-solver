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
  auto working_set = rs.working_set();
  for (auto root : rs.working_set()) {
    candidates[root] = qs.modification_candidates(root);
  }
  return candidates;
}

auto generate_pairs(const std::map<size_t, std::vector<size_t>>& candidates) {
  auto pairs = std::vector<std::pair<size_t, size_t>>{};
  for (const auto& root_candidates : candidates) {
    const auto& c = root_candidates.second;
    for (size_t i = 1; i < c.size(); ++i)
      pairs.emplace_back(i - 1, i);
  }
  return pairs;
};

} // namespace detail

template <class R, class P, class Q, class ST>
void check_conditioning(XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>, ST>& xs, RSpace<R, Q, P>& rs,
                        QSpace<R, Q, P>& qs, PSpace<R, P>& ps, double threshold) {
  bool stable = false;
  auto candidates = detail::generate_candidates(rs, qs);
  while (!stable && !candidates.empty()) {
    auto& s = xs.data[EqnData::S];
    auto svd = svd_system(xs.size(), array::Span<double>{&s(0, 0), s.size()}, threshold);
    stable = svd.empty();
    if (!svd.empty()) {
      auto pairs = detail::generate_pairs(candidates);
      auto norms = std::vector<double>{};
      for (const auto& p : pairs) {
        const auto oQ = xs.dimensions().oQ;
        const auto& v = svd.front().v;
        const auto norm = std::abs(v.at(oQ + p.first) * v.at(oQ + p.second));
        norms.emplace_back(norm);
      }
      auto it_max = std::max_element(begin(norms), end(norms));
      auto i_max = std::distance(begin(norms), it_max);
      qs.merge(pairs.at(i_max));
      xs.build_subspace(rs, qs, ps);
      candidates = detail::generate_candidates(rs, qs);
    }
  }
}

} // namespace xspace
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
