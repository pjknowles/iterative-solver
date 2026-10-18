#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/subspace/DSpace.h>
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>

#include <cassert>
#include <memory>

namespace molpro::linalg::itsolv::subspace {
namespace xspace {
//! New sections of equation data
template <typename T>
struct NewData {
  NewData(size_t nQnew, size_t nX, size_t nRHS) {
    for (auto d : {EqnData::H, EqnData::S}) {
      qq[d].resize({nQnew, nQnew});
      qx[d].resize({nQnew, nX});
      xq[d].resize({nX, nQnew});
    }
    qq[EqnData::rhs].resize({nQnew, nRHS});
  }

  SubspaceData<T> qq = null_data<T, EqnData::H, EqnData::S>(); //!< data block between new paramters
  SubspaceData<T> qx = null_data<T, EqnData::H, EqnData::S>(); //!< data block between new parameters and current X space
  SubspaceData<T> xq = null_data<T, EqnData::H, EqnData::S>(); //!< data block between current X space and new parameters
};

//! Returns new sections of equation data
template <class R, class Q, class P>
auto update_qspace_data(const CVecRef<R>& params, const CVecRef<R>& actions, const CVecRef<P>& pparams,
                        const CVecRef<Q>& qparams, const CVecRef<Q>& qactions, const CVecRef<Q>& dparams,
                        const CVecRef<Q>& dactions, const CVecRef<Q>& rhs, const Dimensions& dims,
                        ArrayHandlers<R, Q, P>& handlers, Logger& logger, bool hermitian = false,
                        bool action_dot_action = false) {
  auto nQnew = params.size();
  auto data = NewData<typename array::ArrayHandler<R, R>::value_type>(nQnew, dims.nX(), rhs.size());
  auto& qq = data.qq;
  auto& qx = data.qx;
  auto& xq = data.xq;
  qq[EqnData::S] = util::overlap(params, handlers.rr());
  qx[EqnData::S].slice({0, dims.oP()}, {nQnew, dims.oP() + dims.nP}) = util::overlap(params, pparams, handlers.rp());
  qx[EqnData::S].slice({0, dims.oQ()}, {nQnew, dims.oQ() + dims.nQ}) = util::overlap(params, qparams, handlers.rq());
  qx[EqnData::S].slice({0, dims.oD()}, {nQnew, dims.oD() + dims.nD}) = util::overlap(params, dparams, handlers.rq());
  qq[EqnData::H] =
      action_dot_action ? util::overlap(actions, handlers.rr()) : util::overlap(params, actions, handlers.rr());
  // left-hand vectors of the H matrix elements with existing actions
  const auto& h_left = action_dot_action ? actions : params;
  qx[EqnData::H].slice({0, dims.oQ()}, {nQnew, dims.oQ() + dims.nQ}) = util::overlap(h_left, qactions, handlers.rq());
  qx[EqnData::H].slice({0, dims.oD()}, {nQnew, dims.oD() + dims.nD}) = util::overlap(h_left, dactions, handlers.rq());
  if (hermitian) {
    xq[EqnData::H].slice({dims.oP(), 0}, {dims.oP() + dims.nP, nQnew}) = util::overlap(pparams, actions, handlers.rp());
    //    xq[EqnData::H].slice({dims.oQ(), 0}, {dims.oQ() + dims.nQ, nQnew}) = util::overlap(action_dot_action ? qactions :
    //    qparams, actions, handlers.qr());
    conjugate_transpose_copy(xq[EqnData::H].slice({dims.oQ(), 0}, {dims.oQ() + dims.nQ, nQnew}),
                   qx[EqnData::H].slice({0, dims.oQ()}, {nQnew, dims.oQ() + dims.nQ}));
    conjugate_transpose_copy(xq[EqnData::H].slice({dims.oD(), 0}, {dims.oD() + dims.nD, nQnew}),
                   qx[EqnData::H].slice({0, dims.oD()}, {nQnew, dims.oD() + dims.nD}));
    conjugate_transpose_copy(qx[EqnData::H].slice({0, dims.oP()}, {nQnew, dims.oP() + dims.nP}),
                   xq[EqnData::H].slice({dims.oP(), 0}, {dims.oP() + dims.nP, nQnew}));
  } else {
    xq[EqnData::H].slice({dims.oQ(), 0}, {dims.oQ() + dims.nQ, nQnew}) = util::overlap(qparams, actions, handlers.rq());
    xq[EqnData::H].slice({dims.oD(), 0}, {dims.oD() + dims.nD, nQnew}) = util::overlap(dparams, actions, handlers.rq());
  }
  qq[EqnData::rhs] = util::overlap(params, rhs, handlers.rq());
  conjugate_transpose_copy(xq[EqnData::S].slice({dims.oP(), 0}, {dims.oP() + dims.nP, nQnew}),
                 qx[EqnData::S].slice({0, dims.oP()}, {nQnew, dims.oP() + dims.nP}));
  conjugate_transpose_copy(xq[EqnData::S].slice({dims.oQ(), 0}, {dims.oQ() + dims.nQ, nQnew}),
                 qx[EqnData::S].slice({0, dims.oQ()}, {nQnew, dims.oQ() + dims.nQ}));
  conjugate_transpose_copy(xq[EqnData::S].slice({dims.oD(), 0}, {dims.oD() + dims.nD, nQnew}),
                 qx[EqnData::S].slice({0, dims.oD()}, {nQnew, dims.oD() + dims.nD}));
  logger.data_dump("xspace::update_qspace_data() nQnew = ", nQnew);
  logger.data_dump("Sqq = ", qq[EqnData::S]);
  logger.data_dump("Hqq = ", qq[EqnData::H]);
  logger.data_dump("Sqx = ", qx[EqnData::S]);
  logger.data_dump("Hqx = ", qx[EqnData::H]);
  logger.data_dump("Sxq = ", xq[EqnData::S]);
  logger.data_dump("Hxq = ", xq[EqnData::H]);
  logger.data_dump("rhs_q = ", qq[EqnData::rhs]);
  return data;
}

//! Calculates overlap blocks between D space and the rest of the subspace
template <class Q, class P>
auto update_dspace_overlap_data(const CVecRef<P>& pparams, const CVecRef<Q>& qparams, const CVecRef<Q>& dparams,
                                const CVecRef<Q>& rhs, array::ArrayHandler<Q, P>& handler_qp,
                                array::ArrayHandler<Q, Q>& handler_qq, Logger& logger) {
  const auto nP = pparams.size();
  const auto nQ = qparams.size();
  const auto nX = nP + nQ;
  auto nD = dparams.size();
  const auto nRHS = rhs.size();
  auto data = NewData<typename array::ArrayHandler<Q, Q>::value_type>(nD, nX, nRHS);
  data.qq[EqnData::S].slice() = util::overlap(dparams, handler_qq);
  data.qx[EqnData::S].slice({0, 0}, {nD, nP}) = util::overlap(dparams, pparams, handler_qp);
  data.qx[EqnData::S].slice({0, nP}, {nD, nX}) = util::overlap(dparams, qparams, handler_qq);
  data.qq[EqnData::rhs].slice() = util::overlap(dparams, rhs, handler_qq);
  conjugate_transpose_copy(data.xq[EqnData::S].slice(), data.qx[EqnData::S].slice());
  logger.data_dump("xspace::update_dspace_overlap_data() nD = ", nD);
  logger.data_dump("Sdd = ", data.qq[EqnData::S]);
  logger.data_dump("Sdx = ", data.qx[EqnData::S]);
  logger.data_dump("rhs_d = ", data.qq[EqnData::rhs]);
  return data;
}

//! Calculates overlap blocks between D space and the rest of the subspace
template <class Q, class P>
auto update_dspace_action_data(const CVecRef<P>& pparams, const CVecRef<Q>& qparams, const CVecRef<Q>& qactions,
                               const CVecRef<Q>& dparams, const CVecRef<Q>& dactions,
                               array::ArrayHandler<Q, P>& handler_qp, array::ArrayHandler<Q, Q>& handler_qq,
                               Logger& logger) {
  const auto nP = pparams.size();
  const auto nQ = qparams.size();
  const auto nX = nP + nQ;
  auto nD = dparams.size();
  auto data = NewData<typename array::ArrayHandler<Q, Q>::value_type>(nD, nX, 0);
  const auto e = EqnData::H;
  data.qq[e].slice() = util::overlap(dparams, dactions, handler_qq);
  data.xq[e].slice({0, 0}, {nP, nD}) = util::overlap(pparams, dactions, handler_qp);
  data.xq[e].slice({nP, 0}, {nX, nD}) = util::overlap(qparams, dactions, handler_qq);
  data.qx[e].slice({0, nP}, {nD, nX}) = util::overlap(dparams, qactions, handler_qq);
  conjugate_transpose_copy(data.qx[e].slice({0, 0}, {nD, nP}), data.xq[e].slice({0, 0}, {nP, nD}));
  logger.data_dump("xspace::update_dspace_action_data() nD = ", nD);
  logger.data_dump("Hdd = ", data.qq[e]);
  logger.data_dump("Hdx = ", data.qx[e]);
  logger.data_dump("Hxd = ", data.xq[e]);
  return data;
}

template <typename T>
void copy_dspace_eqn_data(const NewData<T>& new_data, SubspaceData<T>& data, const subspace::EqnData e,
                          const Dimensions& dims) {
  const auto& dd = new_data.qq.at(e);
  const auto& dx = new_data.qx.at(e);
  const auto& xd = new_data.xq.at(e);
  data[e].slice({dims.oD(), dims.oD()}, {dims.oD() + dims.nD, dims.oD() + dims.nD}) = dd;
  data[e].slice({dims.oD(), dims.oP()}, {dims.oD() + dims.nD, dims.oP() + dims.nP}) =
      dx.slice({0, 0}, {dims.nD, dims.nP});
  data[e].slice({dims.oD(), dims.oQ()}, {dims.oD() + dims.nD, dims.oQ() + dims.nQ}) =
      dx.slice({0, dims.nP}, {dims.nD, dims.nP + dims.nQ});
  data[e].slice({dims.oP(), dims.oD()}, {dims.oP() + dims.nP, dims.oD() + dims.nD}) =
      xd.slice({0, 0}, {dims.nP, dims.nD});
  data[e].slice({dims.oQ(), dims.oD()}, {dims.oQ() + dims.nQ, dims.oD() + dims.nD}) =
      xd.slice({dims.nP, 0}, {dims.nP + dims.nQ, dims.nD});
}
} // namespace xspace

template <class R, class Q, class P>
class XSpace : public IXSpace<R, Q, P> {
public:
  using typename IXSpace<R, Q, P>::value_type;
  using typename IXSpace<R, Q, P>::value_type_abs;
  using IXSpace<R, Q, P>::data;

  explicit XSpace(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers, const std::shared_ptr<Logger>& logger)
      : pspace(), qspace(handlers, logger), dspace(logger), m_handlers(handlers), m_logger(logger) {
    data = null_data<value_type, EqnData::H, EqnData::S, EqnData::rhs>();
  };

  //! Update parameters in Q space and corresponding equation data
  void update_qspace(const CVecRef<R>& params, const CVecRef<R>& actions) override {
    m_logger->trace("QSpace::update_qspace");
    auto new_data =
        xspace::update_qspace_data(params, actions, cparamsp(), cparamsq(), cactionsq(), cparamsd(), cactionsd(),
                                   cwrap(m_rhs), m_dim, *m_handlers, *m_logger, m_hermitian, m_action_dot_action);
    qspace.update(params, actions, new_data.qq, new_data.qx, new_data.xq, m_dim, data);
    update_dimensions();
  }

  //! Clears old D space container and stores new params and actions. @param lin_trans_only_R R space component of D
  void update_dspace(VecRef<Q>& params, VecRef<Q>& actions) override {
    dspace.update(params, actions);
    update_dimensions();
    for (auto e : {EqnData::H, EqnData::S})
      data[e].resize({m_dim.nX(), m_dim.nX()});
    auto new_data = xspace::update_dspace_overlap_data(cparamsp(), cparamsq(), cparamsd(), cwrap(m_rhs),
                                                       m_handlers->qp(), m_handlers->qq(), *m_logger);
    xspace::copy_dspace_eqn_data(new_data, data, EqnData::S, m_dim);
    auto new_data_action = xspace::update_dspace_action_data(
        cparamsp(), cparamsq(), cactionsq(), cparamsd(), cactionsd(), m_handlers->qp(), m_handlers->qq(), *m_logger);
    xspace::copy_dspace_eqn_data(new_data_action, data, EqnData::H, m_dim);
    data[EqnData::rhs].resize({m_dim.nX(), m_dim.nRHS});
    data[EqnData::rhs].slice({m_dim.oD(), 0}, {m_dim.oD() + m_dim.nD, m_dim.nRHS}) = new_data.qq[EqnData::rhs].slice();
  }

  // FIXME this must be called when XSpace is empty
  //! @warning Subspace should be empty
  void update_pspace(const CVecRef<P>& params, const array::Span<value_type>& pp_action_matrix) override {
    assert(m_dim.nX() == 0);
    if (!m_hermitian)
      throw std::runtime_error("P space can only be used with hermitian kernels");
    pspace.update(params, m_handlers->pp());
    update_dimensions();
    const size_t nP = m_dim.nP;
    update_rhs_with_pspace();
    data[EqnData::S].resize({nP, nP});
    data[EqnData::H].resize({nP, nP});
    data[EqnData::S].slice() = util::overlap(params, m_handlers->pp());
    for (size_t i = 0, ij = 0; i < nP; ++i)
      for (size_t j = 0; j < nP; ++j, ++ij)
        data[EqnData::H](i, j) = pp_action_matrix[ij];
  }

  //! For a system of linear equations Ax=b, adds rhs vectors b.
  void add_rhs_equations(const CVecRef<R>& rhs) {
    m_rhs.reserve(m_rhs.size() + rhs.size());
    for (const auto& r : rhs)
      m_rhs.emplace_back(this->m_handlers->qr().copy(r));
    for (const auto& r : rhs) {
      auto d = std::abs(this->m_handlers->rr().dot(r, r));
      if (d == 0)
        throw std::runtime_error("RHS vector cannot be zero");
      m_rhs_norm.emplace_back(std::sqrt(d));
    }
    update_dimensions();
    update_rhs_with_pspace();
  }

  //! Access RHS vectors in linear equations
  CVecRef<Q> rhs() const { return cwrap(m_rhs); }

  //! Norm of RHS vectors
  const std::vector<value_type_abs>& rhs_norm() const { return m_rhs_norm; }

  const Dimensions& dimensions() const override { return m_dim; }

  void erase(size_t i) override {
    if (i >= m_dim.oP() && i < m_dim.oP() + m_dim.nP) {
      erasep(i - m_dim.oP());
    } else if (i >= m_dim.oQ() && i < m_dim.oQ() + m_dim.nQ) {
      eraseq(i - m_dim.oQ());
    } else if (i >= m_dim.oD() && i < m_dim.oD() + m_dim.nD) {
      erased(i - m_dim.oD());
    }
  }

  void eraseq(size_t i) override {
    qspace.erase(i);
    remove_data(m_dim.oQ() + i);
    update_dimensions();
  }

  void erasep(size_t i) override {
    pspace.erase(i);
    remove_data(m_dim.oP() + i);
    update_dimensions();
  }

  void erased(size_t i) override {
    dspace.erase(i);
    remove_data(m_dim.oD() + i);
    update_dimensions();
  }

  VecRef<P> paramsp() override { return pspace.params(); }
  VecRef<Q> paramsq() override { return qspace.params(); }
  VecRef<Q> actionsq() override { return qspace.actions(); }
  VecRef<Q> paramsd() override { return dspace.params(); }
  VecRef<Q> actionsd() override { return dspace.actions(); }

  CVecRef<P> paramsp() const override { return pspace.params(); }
  CVecRef<Q> paramsq() const override { return qspace.params(); }
  CVecRef<Q> actionsq() const override { return qspace.actions(); }
  CVecRef<Q> paramsd() const override { return dspace.cparams(); }
  CVecRef<Q> actionsd() const override { return dspace.cactions(); }

  CVecRef<P> cparamsp() const override { return pspace.cparams(); }
  CVecRef<Q> cparamsq() const override { return qspace.cparams(); }
  CVecRef<Q> cactionsq() const override { return qspace.cactions(); }
  CVecRef<Q> cparamsd() const override { return dspace.cparams(); }
  CVecRef<Q> cactionsd() const override { return dspace.cactions(); }

  void set_logger(std::shared_ptr<Logger> logger) override { m_logger = std::move(logger); }

  //! Set Hermiticity of the subspace. P space can only be used with Hermitian problems
  void set_hermiticity(bool hermitian) { m_hermitian = hermitian; }
  bool get_hermiticity() { return m_hermitian; }
  void set_action_action() { m_action_dot_action = true; }

  PSpace<R, P> pspace;
  QSpace<R, Q, P> qspace;
  DSpace<Q> dspace;

protected:
  void update_dimensions() {
    m_dim = Dimensions(pspace.size(), qspace.size(), dspace.size());
    m_dim.nRHS = m_rhs.size();
  }

  //! Update projection of RHS data onto P space. @warning Subspace should contain only P space
  auto update_rhs_with_pspace() {
    data[EqnData::rhs].resize({m_dim.nP, m_dim.nRHS});
    data[EqnData::rhs].slice() = util::overlap(cparamsp(), rhs(), m_handlers->qp());
  };

  void remove_data(size_t i) {
    for (auto d : {EqnData::H, EqnData::S})
      data[d].remove_row_col(i, i);
    if (data.find(EqnData::rhs) != std::end(data) && !data[EqnData::rhs].empty())
      data[EqnData::rhs].remove_row(i);
    if (data.find(EqnData::value) != std::end(data) && !data[EqnData::value].empty())
      data[EqnData::value].remove_row(i);
  }
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<Logger> m_logger;
  Dimensions m_dim;
  std::vector<Q> m_rhs;                   //!< Right hand side vectors
  std::vector<value_type_abs> m_rhs_norm; //!< norm of RHS vectors
  bool m_hermitian = false;               //!< whether the matrix is Hermitian
  bool m_action_dot_action = false;       //!< whether the H matrix is action.action instead of action.parameter
};

} // namespace molpro::linalg::itsolv::subspace
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
