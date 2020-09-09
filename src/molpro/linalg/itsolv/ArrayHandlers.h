#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_ARRAYHANDLERS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_ARRAYHANDLERS_H
#include <molpro/linalg/array/default_handler.h>

namespace molpro {
namespace linalg {

namespace iterativesolver {
namespace util {

struct ArrayHandlersError : public std::logic_error {
  using std::logic_error::logic_error;
};

} // namespace util

/*!
 * @brief Class, containing a collection of array handlers used in IterativeSolver
 * Provides a Builder sub-class, containing member functions to accumulate user-defined handlers,
 * as well as determine and instantiate default handlers, based on container (template parameter) types
 * @tparam R array for R space
 * @tparam Q array for Q space
 * @tparam P array for P space
 */
template <typename R, typename Q = R, typename P = std::map<size_t, typename R::value_type>>
class ArrayHandlers {
protected:
  class Builder;

public:
  ArrayHandlers(std::shared_ptr<array::ArrayHandler<R, R>> rr, std::shared_ptr<array::ArrayHandler<Q, Q>> qq,
                std::shared_ptr<array::ArrayHandler<P, P>> pp, std::shared_ptr<array::ArrayHandler<R, Q>> rq,
                std::shared_ptr<array::ArrayHandler<R, P>> rp, std::shared_ptr<array::ArrayHandler<Q, R>> qr,
                std::shared_ptr<array::ArrayHandler<Q, P>> qp)
      : m_rr{rr}, m_qq{qq}, m_pp{pp}, m_rq{rq}, m_rp{rp}, m_qr{qr}, m_qp{qp} {}

  //! Uses default handlers
  ArrayHandlers() : ArrayHandlers(Builder{}.build()){};

  auto& rr() { return *m_rr; }
  auto& qq() { return *m_qq; }
  auto& pp() { return *m_pp; }
  auto& rq() { return *m_rq; }
  auto& qr() { return *m_qr; } // TODO qr is not needed, by copy is implemented the wrong way around
  auto& rp() { return *m_rp; }
  auto& qp() { return *m_qp; }

  /*!
   * @brief Utility for creating Array handlers with some user specified handlers
   *
   * For example:
   * @code{c++}
   * // create ArrayHandlers with user defined handlers: rr_explicit, rq_explicit
   * auto array_handlers = ArrayHandlers<Rtype, Qtype, Ptype>::create().rr(rr_explicit).rq(rq_explicit).build();
   * @endcode
   * @return
   */
  static Builder create() { return {}; }

protected:
  class Builder {
  public:
    Builder() : rr(this), qq(this), pp(this), rq(this), rp(this), qr(this), qp(this) {}

    ArrayHandlers build() {
      return ArrayHandlers{rr.handler(), qq.handler(), pp.handler(), rq.handler(),
                           rp.handler(), qr.handler(), qp.handler()};
    };

  protected:
    template <class T, class S>
    class Proxy {
    public:
      Proxy() = default;
      explicit Proxy(Builder* b) : builder(b) {}

      //! assigns a handler to the proxy
      Builder& operator()(const std::shared_ptr<array::ArrayHandler<T, S>>& h) {
        m_handler = h;
        return *builder;
      }

      //! @returns a stored or a default handler
      std::shared_ptr<array::ArrayHandler<T, S>> handler() const {
        if (!m_handler)
          return array::create_default_handler<T, S>();
        else
          return m_handler;
      }

    protected:
      Builder* builder = nullptr;
      std::shared_ptr<array::ArrayHandler<T, S>> m_handler;
    };

  public:
    Proxy<R, R> rr;
    Proxy<Q, Q> qq;
    Proxy<P, P> pp;
    Proxy<R, Q> rq;
    Proxy<R, P> rp;
    Proxy<Q, R> qr;
    Proxy<Q, P> qp;
  };

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
