#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>
#include <molpro/linalg/itsolv/wrap.h>

#include <vector>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace cspace {

//! Update equation data
template <typename T>
void update_subspace(const unsigned int root, SubspaceData& cc, SubspaceData& qc, SubspaceData& cq, SubspaceData& cp,
                     SubspaceData& pc, const T eigval, const SubspaceData& xdata, const Matrix<T>& solutions,
                     const xspace::Dimensions& dims) {
  if (root >= cc[EqnData::S].rows()) {
    for (auto d : {EqnData::S, EqnData::H}) {
      cc[d].resize({root + 1, root + 1});
      qc[d].resize({qc[d].rows(), root + 1});
      cq[d].resize({root + 1, cq[d].cols()});
      pc[d].resize({qc[d].rows(), root + 1});
      cp[d].resize({root + 1, cp[d].cols()});
    }
  }
  cc[EqnData::S](root, root) = 1;
  cc[EqnData::H](root, root) = eigval;
  auto update_offdiag = [&solutions, &dims, &xdata, root](SubspaceData& cy, SubspaceData& yc, size_t oY, size_t nY) {
    for (size_t i = 0; i < nY; ++i) {
      cy[EqnData::S](root, i) = yc[EqnData::S](i, root) = solutions(root, oY + i);
      cy[EqnData::H](root, i) = yc[EqnData::H](i, root) = 0;
      for (size_t j = 0; j < dims.nX; ++j) {
        cy[EqnData::H](root, i) += solutions(root, j) * xdata.at(EqnData::H)(j, oY + i);
        yc[EqnData::H](i, root) += xdata.at(EqnData::H)(oY + i, j) * solutions(root, j);
      }
    }
  };
  update_offdiag(cp, pc, dims.oP, dims.nP);
  update_offdiag(cq, qc, dims.oQ, dims.nQ);
}

} // namespace cspace

/*!
 * @brief Space storing the complement of Q necessary to reconstruct previous solutions.
 *
 * It consists of linear combinations of deleted Q parameters and maximum size of the space is the number of solutions.
 * As such, C space accumulates round off error and eventually actions will not map to corresponding parameters
 * leading to errors in the subspace. Thus C space should be periodically reset.
 *
 * Resetting C space
 * -----------------
 *  Use C space as new parameters, calculate their exact action and add to Q.
 *
 *  Cons: size of C space can be greater than R, so what will happen to the rest of C space?
 *
 * This can be done over multiple iterations, as long as no more Q parameters are removed. Thus firstly remove nC
 * parameters from Q, than start moving C parameters to R params in propose_rspace, until C space is empty. At that
 * point all of C space was moved into Q and there were no modifications, since size of Q space was within limit and
 * residual was not used. It should be possible to have a mixed update where new parameters are orthogonalised residuals
 * and C params.
 */
template <class Rt, class Qt, class Pt, typename ST>
class CSpace {
public:
  using R = Rt;
  using Q = Qt;
  using P = Pt;
  using scalar_type = ST;

  //! Matrix and overlap data mapped to the subspace
  SubspaceData data = null_data<EqnData::H, EqnData::S>();
  SubspaceData qc = null_data<EqnData::H, EqnData::S>();
  SubspaceData cq = null_data<EqnData::H, EqnData::S>();
  SubspaceData cp = null_data<EqnData::H, EqnData::S>();
  SubspaceData pc = null_data<EqnData::H, EqnData::S>();

  explicit CSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_logger(std::move(logger)) {}

  //! Adds a solution to the space, overwriting if a solution for the root already exists, and construct equation data
  void update(const unsigned int root, const R& param, const R& action, scalar_type error) {
    auto it = m_params.find(root);
    if (it == m_params.end()) {
      m_params.emplace(root, m_handlers->qr().copy(param));
      m_actions.emplace(root, m_handlers->qr().copy(action));
      m_errors[root] = error;
    } else {
      m_handlers->qr().copy(m_params.at(root), param);
      m_handlers->qr().copy(m_actions.at(root), action);
      m_errors[root] = error;
    }
  }

  //! Adds solutions to the space, overwriting if a solution for given root already exists, and construct equation data
  void update(const std::vector<unsigned int>& roots, const CVecRef<R>& params, const CVecRef<R>& actions,
              const std::vector<scalar_type>& errors) {
    for (size_t i = 0; i < roots.size(); ++i) {
      update(roots.at(i), params.at(i), actions.at(i), errors.at(i));
    }
  }

  void set_error(unsigned int root, scalar_type error) { m_errors[root] = error; }
  void set_error(const std::vector<unsigned int>& roots, const std::vector<scalar_type>& errors) {
    for (size_t i = 0; i < roots.size(); ++i) {
      set_error(roots[i], errors[i]);
    }
  }

  size_t size() const { return m_params.size(); }

  //! Removes all vectors and data
  void clear() {
    m_params.clear();
    m_actions.clear();
    m_errors.clear();
    for (auto d : data) {
      d.second.clear();
    }
  }

  //! Erases parameter i. @param i index in the current space
  void erase(size_t i) {
    assert(m_params.size() > i);
    auto erase_at_i = [i](auto& v) { v.erase(std::next(begin(v), i)); };
    erase_at_i(m_params);
    erase_at_i(m_actions);
    erase_at_i(m_errors);
    for (auto d : {EqnData::H, EqnData::S}) {
      data.at(d).remove_row_col(i, i);
    }
  }

  //! Number of solutions stored in C space
  size_t nroots() const {
    if (m_errors.empty())
      return 0;
    else
      return m_errors.rbegin().base()->first;
  }

  CVecRef<Q> params() const { return wrap(m_params); };
  CVecRef<Q> cparams() const { return params(); };
  VecRef<Q> params() { return wrap(m_params); };
  CVecRef<Q> actions() const { return wrap(m_actions); };
  CVecRef<Q> cactions() const { return actions(); };
  VecRef<Q> actions() { return wrap(m_actions); };

  std::vector<scalar_type> errors() const {
    auto err = std::vector<scalar_type>(nroots());
    for (const auto& e : m_errors)
      err.at(e.first) = e.second;
    return err;
  };

protected:
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<Logger> m_logger;
  std::map<unsigned int, Q> m_params;           //! parameters for root index
  std::map<unsigned int, Q> m_actions;          //! actions for root index
  std::map<unsigned int, scalar_type> m_errors; //! errors for root index
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
