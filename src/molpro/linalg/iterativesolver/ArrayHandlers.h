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
public:
  ArrayHandlers(std::shared_ptr<array::ArrayHandler<R, R>> rr, std::shared_ptr<array::ArrayHandler<Q, Q>> qq,
                std::shared_ptr<array::ArrayHandler<P, P>> pp, std::shared_ptr<array::ArrayHandler<R, Q>> rq,
                std::shared_ptr<array::ArrayHandler<R, P>> rp, std::shared_ptr<array::ArrayHandler<Q, R>> qr,
                std::shared_ptr<array::ArrayHandler<Q, P>> qp)
      : m_rr{rr}, m_qq{qq}, m_pp{pp}, m_rq{rq}, m_rp{rp}, m_qr{qr}, m_qp{qp} {}

  class Builder {

  public:
    Builder& rr(std::shared_ptr<array::ArrayHandler<R, R>> handler) {
      m_rr = handler;
      return *this;
    }

    Builder& qq(std::shared_ptr<array::ArrayHandler<Q, Q>> handler) {
      m_qq = handler;
      return *this;
    }

    Builder& pp(std::shared_ptr<array::ArrayHandler<P, P>> handler) {
      m_pp = handler;
      return *this;
    }

    Builder& rq(std::shared_ptr<array::ArrayHandler<R, Q>> handler) {
      m_rq = handler;
      return *this;
    }

    Builder& rp(std::shared_ptr<array::ArrayHandler<R, P>> handler) {
      m_rp = handler;
      return *this;
    }

    Builder& qr(std::shared_ptr<array::ArrayHandler<Q, R>> handler) {
      m_qr = handler;
      return *this;
    }

    Builder& qp(std::shared_ptr<array::ArrayHandler<Q, P>> handler) {
      m_qp = handler;
      return *this;
    }

    void add_default_handlers() const {
      if (!m_rr) {
        m_rr = array::create_default_handler<R, R>();
      }
      if (!m_qq) {
        m_qq = array::create_default_handler<Q, Q>();
      }
      if (!m_pp) {
        m_pp = array::create_default_handler<P, P>();
      }
      if (!m_rq) {
        m_rq = array::create_default_handler<R, Q>();
      }
      if (!m_rp) {
        m_rp = array::create_default_handler<R, P>();
      }
      if (!m_qr) {
        m_qr = array::create_default_handler<Q, R>();
      }
      if (!m_qp) {
        m_qp = array::create_default_handler<Q, P>();
      }
    }

    ArrayHandlers build() const {
      add_default_handlers();
      return ArrayHandlers{m_rr, m_qq, m_pp, m_rq, m_rp, m_qr, m_qp};
    };

  private:
    mutable std::shared_ptr<array::ArrayHandler<R, R>> m_rr;
    mutable std::shared_ptr<array::ArrayHandler<Q, Q>> m_qq;
    mutable std::shared_ptr<array::ArrayHandler<P, P>> m_pp;
    mutable std::shared_ptr<array::ArrayHandler<R, Q>> m_rq;
    mutable std::shared_ptr<array::ArrayHandler<R, P>> m_rp;
    mutable std::shared_ptr<array::ArrayHandler<Q, R>> m_qr;
    mutable std::shared_ptr<array::ArrayHandler<Q, P>> m_qp;

    void error(const std::string& message) { throw util::ArrayHandlersError{message}; };
  };

  auto& rr() { return *m_rr; }
  auto& qq() { return *m_qq; }
  auto& pp() { return *m_pp; }
  auto& rq() { return *m_rq; }
  auto& qr() { return *m_qr; } // TODO qr is not needed, by copy is implemented the wrong way around
  auto& rp() { return *m_rp; }
  auto& qp() { return *m_qp; }

private:
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
