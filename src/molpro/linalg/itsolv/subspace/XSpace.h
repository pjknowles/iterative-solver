#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
#include <cassert>
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/DSpace.h>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
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
                        const CVecRef<Q>& qparams, const CVecRef<Q>& qactions, const CVecRef<Q>& dparams,
                        const CVecRef<Q>& dactions, const Dimensions& dims, ArrayHandlers<R, Q, P>& handlers,
                        Logger& logger) {
  auto nQnew = params.size();
  auto data = NewData(nQnew, dims.nX);
  auto& qq = data.qq;
  auto& qx = data.qx;
  auto& xq = data.xq;
  qq[EqnData::S] = util::overlap(params, handlers.rr());
  qx[EqnData::S].slice({0, dims.oP}, {nQnew, dims.oP + dims.nP}) = util::overlap(params, pparams, handlers.rp());
  qx[EqnData::S].slice({0, dims.oQ}, {nQnew, dims.oQ + dims.nQ}) = util::overlap(params, qparams, handlers.rq());
  qx[EqnData::S].slice({0, dims.oD}, {nQnew, dims.oD + dims.nD}) = util::overlap(params, dparams, handlers.rq());
  qq[EqnData::H] = util::overlap(params, actions, handlers.rr());
  qx[EqnData::H].slice({0, dims.oQ}, {nQnew, dims.oQ + dims.nQ}) = util::overlap(params, qactions, handlers.rq());
  qx[EqnData::H].slice({0, dims.oD}, {nQnew, dims.oD + dims.nD}) = util::overlap(params, dactions, handlers.rq());
  xq[EqnData::H].slice({dims.oP, 0}, {dims.oP + dims.nP, nQnew}) = util::overlap(pparams, actions, handlers.rp());
  xq[EqnData::H].slice({dims.oQ, 0}, {dims.oQ + dims.nQ, nQnew}) = util::overlap(qparams, actions, handlers.qr());
  xq[EqnData::H].slice({dims.oD, 0}, {dims.oD + dims.nD, nQnew}) = util::overlap(dparams, actions, handlers.qr());
  transpose_copy(xq[EqnData::S].slice({dims.oP, 0}, {dims.oP + dims.nP, nQnew}),
                 qx[EqnData::S].slice({0, dims.oP}, {nQnew, dims.oP + dims.nP}));
  transpose_copy(xq[EqnData::S].slice({dims.oQ, 0}, {dims.oQ + dims.nQ, nQnew}),
                 qx[EqnData::S].slice({0, dims.oQ}, {nQnew, dims.oQ + dims.nQ}));
  transpose_copy(xq[EqnData::S].slice({dims.oD, 0}, {dims.oD + dims.nD, nQnew}),
                 qx[EqnData::S].slice({0, dims.oD}, {nQnew, dims.oD + dims.nD}));
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
} // namespace xspace

template <class R, class Q, class P>
class XSpace : public XSpaceI<R, Q, P> {
public:
  using typename XSpaceI<R, Q, P>::value_type;
  using typename XSpaceI<R, Q, P>::value_type_abs;
  using XSpaceI<R, Q, P>::data;

  explicit XSpace(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers, const std::shared_ptr<Logger>& logger)
      : pspace(), qspace(handlers, logger), dspace(logger), m_handlers(handlers), m_logger(logger) {
    data = null_data<EqnData::H, EqnData::S>();
  };

  //! Updata parameters in Q space and corresponding equation data
  void update_qspace(const CVecRef<R>& params, const CVecRef<R>& actions) override {
    m_logger->msg("QSpace::update_qspace", Logger::Trace);
    auto new_data = xspace::update_qspace_data(params, actions, cparamsp(), cparamsq(), cactionsq(), cparamsd(),
                                               cactionsd(), m_dim, *m_handlers, *m_logger);
    qspace.update(params, actions, new_data.qq, new_data.qx, new_data.xq, m_dim, data);
    update_dimensions();
  }

  //! Clears old D space container and stores new params and actions. @param lin_trans_only_R R space component of D
  void update_dspace(VecRef<R>& params, VecRef<R>& actions, Matrix<value_type>& lin_trans_only_R) override {
    dspace.update(params, actions, lin_trans_only_R);
    update_dimensions();
    for (auto e : {EqnData::H, EqnData::S})
      data[e].resize({m_dim.nX, m_dim.nX});
  }

  const xspace::Dimensions& dimensions() const override { return m_dim; }

  void erase(size_t i) override {
    if (m_dim.oP >= i && i < m_dim.oP + m_dim.nP) {
      erasep(i - m_dim.oP);
    } else if (m_dim.oQ >= i && i < m_dim.oQ + m_dim.nQ) {
      eraseq(i - m_dim.oQ);
    } else if (m_dim.oD >= i && i < m_dim.oD + m_dim.nD) {
      erased(i - m_dim.oD);
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

  void erased(size_t i) override {
    dspace.erase(i);
    remove_data(m_dim.oD + i);
    update_dimensions();
  }

  VecRef<P> paramsp() override { return pspace.params(); }
  VecRef<P> actionsp() override { return pspace.actions(); }
  VecRef<Q> paramsq() override { return qspace.params(); }
  VecRef<Q> actionsq() override { return qspace.actions(); }
  VecRef<Q> paramsd() override { return dspace.params(); }
  VecRef<Q> actionsd() override { return dspace.actions(); }

  CVecRef<P> paramsp() const override { return pspace.params(); }
  CVecRef<P> actionsp() const override { return pspace.actions(); }
  CVecRef<Q> paramsq() const override { return qspace.params(); }
  CVecRef<Q> actionsq() const override { return qspace.actions(); }
  CVecRef<Q> paramsd() const override { return dspace.cparams(); }
  CVecRef<Q> actionsd() const override { return dspace.cactions(); }

  CVecRef<P> cparamsp() const override { return pspace.cparams(); }
  CVecRef<P> cactionsp() const override { return pspace.cactions(); }
  CVecRef<Q> cparamsq() const override { return qspace.cparams(); }
  CVecRef<Q> cactionsq() const override { return qspace.cactions(); }
  CVecRef<Q> cparamsd() const override { return dspace.cparams(); }
  CVecRef<Q> cactionsd() const override { return dspace.cactions(); }

  PSpace<R, P> pspace;
  QSpace<R, Q, P> qspace;
  DSpace<Q> dspace;

protected:
  void update_dimensions() { m_dim = xspace::Dimensions(pspace.size(), qspace.size(), dspace.size()); }

  void remove_data(size_t i) {
    for (auto d : {EqnData::H, EqnData::S})
      data[d].remove_row_col(i, i);
  }
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<Logger> m_logger;
  xspace::Dimensions m_dim;
  bool m_hermitian = false; //!< whether the matrix is Hermitian
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
