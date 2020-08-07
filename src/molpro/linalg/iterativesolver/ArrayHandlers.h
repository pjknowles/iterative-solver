#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_ARRAYHANDLERS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_ARRAYHANDLERS_H
#include <molpro/linalg/array/ArrayHandlerFactory.h>

namespace molpro {
namespace linalg {
namespace iterativesolver {

// TODO We need a constructor that can take any number of handlers explicitly and generate the rest using the factory
/*!
 * @brief Collection of array handlers needed by IterativeSolver
 * @tparam R array for R space
 * @tparam Q array for Q space
 * @tparam P array for P space
 */
template <typename R, typename Q = R, typename P = std::map<size_t, typename R::value_type>>
struct ArrayHandlers {
  auto& rr() { return *m_rr; }
  auto& qq() { return *m_qq; }
  auto& pp() { return *m_pp; }
  auto& rq() { return *m_rq; }
  auto& qr() { return *m_qr; } // TODO qr is not needed, by copy is implemented the wrong way around
  auto& rp() { return *m_rp; }
  auto& qp() { return *m_qp; }

  std::shared_ptr<array::ArrayHandler<R, R>> m_rr;
  std::shared_ptr<array::ArrayHandler<Q, Q>> m_qq;
  std::shared_ptr<array::ArrayHandler<P, P>> m_pp;
  std::shared_ptr<array::ArrayHandler<R, Q>> m_rq;
  std::shared_ptr<array::ArrayHandler<R, P>> m_rp;
  std::shared_ptr<array::ArrayHandler<Q, R>> m_qr;
  std::shared_ptr<array::ArrayHandler<Q, P>> m_qp;
};
} // namespace iterativesolver
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_ARRAYHANDLERS_H
