#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
#include <cassert>
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/CSpace.h>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>

namespace molpro::linalg::itsolv::subspace {
namespace xspace {
//! New sections of equation data
struct NewData {
  NewData(size_t nQnew, size_t nX) {
    for (auto d : {EqnData::H, EqnData::S}) {
      qq[d].resize({nQnew, nQnew});
      qx[d].resize({nQnew, nX});
      xq[d].resize({nX, nQnew});
    }
  }

  SubspaceData qq = null_data<EqnData::H, EqnData::S>(); //!< data block between new paramters
  SubspaceData qx = null_data<EqnData::H, EqnData::S>(); //!< data block between new parameters and current X space
  SubspaceData xq = null_data<EqnData::H, EqnData::S>(); //!< data block between current X space and new parameters
};

//! Returns new sections of equation data
template <class R, class Q, class P>
auto update_qspace_data(const CVecRef<R>& params, const CVecRef<R>& actions, const CVecRef<P>& pparams,
                        const CVecRef<Q>& qparams, const CVecRef<Q>& qactions, const CVecRef<Q>& cparams,
                        const CVecRef<Q>& cactions, const Dimensions& dims, ArrayHandlers<R, Q, P>& handlers,
                        Logger& logger) {
  auto nQnew = params.size();
  auto data = NewData(nQnew, dims.nX);
  auto& qq = data.qq;
  auto& qx = data.qx;
  auto& xq = data.xq;
  qq[EqnData::S] = util::overlap(params, handlers.rr());
  qx[EqnData::S].slice({0, dims.oP}, {nQnew, dims.oP + dims.nP}) = util::overlap(params, pparams, handlers.rp());
  qx[EqnData::S].slice({0, dims.oQ}, {nQnew, dims.oQ + dims.nQ}) = util::overlap(params, qparams, handlers.rq());
  qx[EqnData::S].slice({0, dims.oC}, {nQnew, dims.oC + dims.nC}) = util::overlap(params, cparams, handlers.rq());
  qq[EqnData::H] = util::overlap(params, actions, handlers.rr());
  qx[EqnData::H].slice({0, dims.oQ}, {nQnew, dims.oQ + dims.nQ}) = util::overlap(params, qactions, handlers.rq());
  qx[EqnData::H].slice({0, dims.oC}, {nQnew, dims.oC + dims.nC}) = util::overlap(params, cactions, handlers.rq());
  xq[EqnData::H].slice({dims.oP, 0}, {dims.oP + dims.nP, nQnew}) = util::overlap(pparams, actions, handlers.rp());
  xq[EqnData::H].slice({dims.oQ, 0}, {dims.oQ + dims.nQ, nQnew}) = util::overlap(qparams, actions, handlers.qr());
  xq[EqnData::H].slice({dims.oC, 0}, {dims.oC + dims.nC, nQnew}) = util::overlap(cparams, actions, handlers.qr());
  transpose_copy(xq[EqnData::S].slice({dims.oP, 0}, {dims.oP + dims.nP, nQnew}),
                 qx[EqnData::S].slice({0, dims.oP}, {nQnew, dims.oP + dims.nP}));
  transpose_copy(xq[EqnData::S].slice({dims.oQ, 0}, {dims.oQ + dims.nQ, nQnew}),
                 qx[EqnData::S].slice({0, dims.oQ}, {nQnew, dims.oQ + dims.nQ}));
  transpose_copy(xq[EqnData::S].slice({dims.oC, 0}, {dims.oC + dims.nC, nQnew}),
                 qx[EqnData::S].slice({0, dims.oC}, {nQnew, dims.oC + dims.nC}));
  // FIXME only works for Hermitian Hamiltonian
  transpose_copy(qx[EqnData::H].slice({0, dims.oP}, {nQnew, dims.oP + dims.nP}),
                 xq[EqnData::H].slice({dims.oP, 0}, {dims.oP + dims.nP, nQnew}));
  if (logger.data_dump) {
    logger.msg("qspace::update_qspace_data() nQnew = " + std::to_string(nQnew), Logger::Info);
    logger.msg("Sqq = " + as_string(qq[EqnData::S]), Logger::Info);
    logger.msg("Hqq = " + as_string(qq[EqnData::H]), Logger::Info);
    logger.msg("Sqx = " + as_string(qx[EqnData::S]), Logger::Info);
    logger.msg("Hqx = " + as_string(qx[EqnData::H]), Logger::Info);
    logger.msg("Sxq = " + as_string(xq[EqnData::S]), Logger::Info);
    logger.msg("Hxq = " + as_string(xq[EqnData::H]), Logger::Info);
  }
  return data;
}

/*!
 * @brief Constructs data blocks for new solutions
 *
 * Equation data is constructed from the subspace data instead of calculating overlaps among full parameters.
 * Both approaches would lead to the same round-off error, but working in the subspace is much cheaper.
 *
 * Solutions {c} are expressed in the current subspace {x}
 * c_i = \sum_j C_{i,j} x_j
 *
 * S_{i,j} = <c_i, x_j> = \sum_k C_{i,k} <x_k, x_j>
 * H_{i,j} = <c_i, x_j>_H = \sum_k C_{i,k} <x_k, x_j>_H
 *
 * <x_k, x_j>_H = <x_k| H |x_j>
 */
template <typename value_type>
auto update_cspace_data(const Matrix<value_type>& solutions, const std::vector<value_type>& eigenvalues,
                        SubspaceData& data, const Dimensions& dims, Logger& logger) {
  const auto nCnew = solutions.rows();
  auto new_data_blocks = NewData{nCnew, dims.nX};
  auto& cc = new_data_blocks.qq;
  auto& cx = new_data_blocks.qx;
  auto& xc = new_data_blocks.xq;
  for (size_t i = 0; i < nCnew; ++i) {
    cc[EqnData::S](i, i) = 1;
    cc[EqnData::H](i, i) = eigenvalues.at(i);
  }
  for (size_t i = 0; i < nCnew; ++i) {
    for (size_t j = 0; j < dims.nX; ++j) {
      for (size_t k = 0; k < dims.nX; ++k) {
        cx[EqnData::S](i, j) += solutions(i, k) * data.at(EqnData::S)(k, j);
        cx[EqnData::H](i, j) += solutions(i, k) * data.at(EqnData::H)(k, j);
        xc[EqnData::H](j, i) += data.at(EqnData::H)(j, k) * solutions(i, k);
      }
      xc[EqnData::S](j, i) = cx[EqnData::S](i, j);
    }
  }
  if (logger.data_dump) {
    logger.msg("qspace::update_cspace_data() nCnew = " + std::to_string(nCnew), Logger::Info);
    logger.msg("Scc = " + as_string(cc[EqnData::S]), Logger::Info);
    logger.msg("Hcc = " + as_string(cc[EqnData::H]), Logger::Info);
    logger.msg("Scx = " + as_string(cx[EqnData::S]), Logger::Info);
    logger.msg("Hcx = " + as_string(cx[EqnData::H]), Logger::Info);
    logger.msg("Sxc = " + as_string(xc[EqnData::S]), Logger::Info);
    logger.msg("Hxc = " + as_string(xc[EqnData::H]), Logger::Info);
  }
  auto new_data = null_data<EqnData::S, EqnData::H>();
  const auto nXnew = dims.nX + nCnew - dims.nC;
  for (auto d : {EqnData::S, EqnData::H}) {
    new_data[d].resize({nXnew, nXnew});
    new_data[d].slice({0, 0}, {dims.oC, dims.oC}) = data.at(d).slice({0, 0}, {dims.oC, dims.oC});
    new_data[d].slice({dims.oC + nCnew, dims.oC + nCnew}, {nXnew, nXnew}) =
        data.at(d).slice({dims.oC + dims.nC, dims.oC + dims.nC}, {dims.nX, dims.nX});
    new_data[d].slice({dims.oC, dims.oC}, {dims.oC + nCnew, dims.oC + nCnew}) = cc.at(d);
    new_data[d].slice({dims.oC, 0}, {dims.oC + nCnew, dims.oC}) = cx.at(d).slice({0, 0}, {nCnew, dims.oC});
    new_data[d].slice({dims.oC, dims.oC + nCnew}, {dims.oC + nCnew, nXnew}) =
        cx.at(d).slice({0, dims.oC + dims.nC}, {nCnew, dims.nX});
    new_data[d].slice({0, dims.oC}, {dims.oC, dims.oC + nCnew}) = xc.at(d).slice({0, 0}, {dims.oC, nCnew});
    new_data[d].slice({dims.oC + nCnew, dims.oC}, {nXnew, dims.oC + nCnew}) =
        xc.at(d).slice({dims.oC + dims.nC, 0}, {dims.nX, nCnew});
    data[d] = std::move(new_data[d]);
  }
}

} // namespace xspace

template <class R, class Q, class P>
class XSpace : public XSpaceI<R, Q, P> {
public:
  using typename XSpaceI<R, Q, P>::value_type;
  using typename XSpaceI<R, Q, P>::value_type_abs;
  using XSpaceI<R, Q, P>::data;

  explicit XSpace(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers, const std::shared_ptr<Logger>& logger)
      : pspace(), qspace(handlers, logger), cspace(handlers, logger), m_handlers(handlers), m_logger(logger) {
    data = null_data<EqnData::H, EqnData::S>();
  };

  //! Updata parameters in Q space and corresponding equation data
  void update_qspace(const CVecRef<R>& params, const CVecRef<R>& actions) override {
    m_logger->msg("QSpace::update_qspace", Logger::Trace);
    auto new_data = xspace::update_qspace_data(params, actions, cparamsp(), cparamsq(), cactionsq(), cparamsc(),
                                               cactionsc(), m_dim, *m_handlers, *m_logger);
    qspace.update(params, actions, new_data.qq, new_data.qx, new_data.xq, m_dim, data);
    update_dimensions();
  }

  //! Uses solutions to update equation data in the subspace
  void update_cspace_data(const Matrix<value_type>& evec, const std::vector<value_type>& eval) {
    m_logger->msg("QSpace::update_cspace_data", Logger::Trace);
    xspace::update_cspace_data(evec, eval, data, m_dim, *m_logger);
  }

  //! Updates a solution in C space. Equation data is not modified (see update_cspace_data())
  void update_cspace(const std::vector<unsigned int>& roots, const CVecRef<R>& params, const CVecRef<R>& actions,
                     const std::vector<value_type>& errors) override {
    m_logger->msg("QSpace::update_cspace", Logger::Trace);
    cspace.update(roots, params, actions, errors);
    update_dimensions();
  }

  const xspace::Dimensions& dimensions() const override { return m_dim; }

  void erase(size_t i) override {
    if (m_dim.oP >= i && i < m_dim.oP + m_dim.nP) {
      erasep(i - m_dim.oP);
    } else if (m_dim.oQ >= i && i < m_dim.oQ + m_dim.nQ) {
      eraseq(i - m_dim.oQ);
    } else if (m_dim.oC >= i && i < m_dim.oC + m_dim.nC) {
      erasec(i - m_dim.oC);
    }
  }

  void eraseq(size_t i) override {
    qspace.erase(i);
    remove_data(m_dim.oQ + i);
    update_dimensions();
  }

  void erasep(size_t i) override {
    pspace.erase(i);
    remove_data(m_dim.oP + i);
    update_dimensions();
  }

  void erasec(size_t i) override {
    cspace.erase(i);
    remove_data(m_dim.oC + i);
    update_dimensions();
  }

  VecRef<P> paramsp() override { return pspace.params(); }
  VecRef<P> actionsp() override { return pspace.actions(); }
  VecRef<Q> paramsq() override { return qspace.params(); }
  VecRef<Q> actionsq() override { return qspace.actions(); }
  VecRef<Q> paramsc() override { return cspace.params(); }
  VecRef<Q> actionsc() override { return cspace.actions(); }

  CVecRef<P> paramsp() const override { return pspace.params(); }
  CVecRef<P> actionsp() const override { return pspace.actions(); }
  CVecRef<Q> paramsq() const override { return qspace.params(); }
  CVecRef<Q> actionsq() const override { return qspace.actions(); }
  CVecRef<Q> paramsc() const override { return cspace.cparams(); }
  CVecRef<Q> actionsc() const override { return cspace.cactions(); }

  CVecRef<P> cparamsp() const override { return pspace.cparams(); }
  CVecRef<P> cactionsp() const override { return pspace.cactions(); }
  CVecRef<Q> cparamsq() const override { return qspace.cparams(); }
  CVecRef<Q> cactionsq() const override { return qspace.cactions(); }
  CVecRef<Q> cparamsc() const override { return cspace.cparams(); }
  CVecRef<Q> cactionsc() const override { return cspace.cactions(); }

  PSpace<R, P> pspace;
  QSpace<R, Q, P> qspace;
  CSpace<R, Q, P, value_type> cspace;

protected:
  void update_dimensions() { m_dim = xspace::Dimensions(pspace.size(), qspace.size(), cspace.size()); }

  void remove_data(size_t i) {
    for (auto d : {EqnData::H, EqnData::S})
      data[d].remove_row_col(i, i);
  }
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<Logger> m_logger;
  xspace::Dimensions m_dim;
  bool m_hermitian = false; //!< whether the matrix is Hermitian
};

} // namespace molpro::linalg::itsolv::subspace
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
