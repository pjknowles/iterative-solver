#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_ARRAYHANDLERS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_ARRAYHANDLERS_H
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/ArrayHandlerSparse.h>
#include <molpro/linalg/array/ArrayHandlerIterable.h>
//#include <molpro/linalg/array/ArrayHandlerIterableSparse.h>
#include <molpro/linalg/array/ArrayHandlerDistr.h>
//#include <molpro/linalg/array/ArrayHandlerDistrSparse.h>

namespace molpro {
namespace linalg {
namespace iterativesolver {
namespace util {

struct ArrayHandlersError : public std::logic_error {
  using std::logic_error::logic_error;
};

}
// TODO We need a constructor that can take any number of handlers explicitly and generate the rest using the factory
/*!
 * @brief Collection of array handlers needed by IterativeSolver
 * @tparam R array for R space
 * @tparam Q array for Q space
 * @tparam P array for P space
 */
//template <typename R, typename Q = R, typename P = std::map<size_t, typename R::value_type>>
//struct ArrayHandlers {
//  auto& rr() { return *m_rr; }
//  auto& qq() { return *m_qq; }
//  auto& pp() { return *m_pp; }
//  auto& rq() { return *m_rq; }
//  auto& qr() { return *m_qr; } // TODO qr is not needed, by copy is implemented the wrong way around
//  auto& rp() { return *m_rp; }
//  auto& qp() { return *m_qp; }
//
//  std::shared_ptr<array::ArrayHandler<R, R>> m_rr;
//  std::shared_ptr<array::ArrayHandler<Q, Q>> m_qq;
//  std::shared_ptr<array::ArrayHandler<P, P>> m_pp;
//  std::shared_ptr<array::ArrayHandler<R, Q>> m_rq;
//  std::shared_ptr<array::ArrayHandler<R, P>> m_rp;
//  std::shared_ptr<array::ArrayHandler<Q, R>> m_qr;
//  std::shared_ptr<array::ArrayHandler<Q, P>> m_qp;
//};

/*!
 * @brief Class, containing a collection of array handlers used in IterativeSolver
 * Provides a Builder sub-class, containing member functions to accumulate user-defined handlers,
 * as well as determine and instantiate default handlers, based on container (template parameter) types
 */
template <typename R, typename Q = R, typename P = std::map<size_t, typename R::value_type>>
class ArrayHandlers {
public:

  ArrayHandlers(std::shared_ptr<array::ArrayHandler<R, R>> rr, std::shared_ptr<array::ArrayHandler<Q, Q>> qq,
                std::shared_ptr<array::ArrayHandler<P, P>> pp, std::shared_ptr<array::ArrayHandler<R, Q>> rq,
                std::shared_ptr<array::ArrayHandler<R, P>> rp, std::shared_ptr<array::ArrayHandler<Q, R>> qr,
                std::shared_ptr<array::ArrayHandler<Q, P>> qp) : m_rr{rr}, m_qq{qq}, m_pp{pp}, m_rq{rq}, m_rp{rp},
                                                                 m_qr{qr}, m_qp{qp} {}

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
    
    //template <typename S, typename T = S>
    //Builder& add_my_handler(std::shared_ptr<array::ArrayHandler<S, T>> handler) {
    //  map_handler<S, T>(handler);
    //  //map_handler(handler);
    //  return *this;
    //};

    //template<typename S, typename T = S>
    //Builder& add_default_handler_iterable() {
    //  std::shared_ptr<array::ArrayHandlerIterable<S, T>> handler = std::make_shared<array::ArrayHandlerIterable<S, T>>();
    //  map_handler<S, T>(handler);
    //  //map_handler(handler);
    //  return *this;
    //};
    
    void add_default_handlers() const {
      if (!m_rr) {
        //m_rr = create_default_handler<R, R>();
        m_rr = std::make_shared<array::ArrayHandlerIterable<R, R>>();
      }
      if (!m_qq) {
        //m_qq = create_default_handler<Q, Q>();
        m_qq = std::make_shared<array::ArrayHandlerIterable<Q, Q>>();
      }
      if (!m_pp) {
        m_pp = create_default_handler<P, P>();
        //m_pp = std::make_shared<array::ArrayHandlerIterable<P, P>>();
      }
      if (!m_rq) {
        //m_rq = create_default_handler<R, Q>();
        m_rq = std::make_shared<array::ArrayHandlerIterable<R, Q>>();
      }
      //if (!m_rp) {
      //  //m_rp = create_default_handler<R, P>();
      //  m_rp = std::make_shared<array::ArrayHandlerIterable<R, P>>();
      //}
      if (!m_qr) {
        //m_qr = create_default_handler<Q, R>();
        m_qr = std::make_shared<array::ArrayHandlerIterable<Q, R>>();
      }
      //if (!m_qp) {
      //  //m_qp = create_default_handler<Q, P>();
      //  m_qp = std::make_shared<array::ArrayHandlerIterable<Q, P>>();
      //}
    }
    
    template< typename S, typename T,
              typename = std::enable_if_t<array::util::has_mapped_type<S>{}>,
              typename = std::enable_if_t<array::util::has_mapped_type<T>{}> >
    static auto create_default_handler() {
      return std::make_shared<array::ArrayHandlerSparse<S,T>>();
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
  
    //template<typename S, typename T = S>
    //void map_handler(std::shared_ptr<array::ArrayHandler<S, T>> handler) {
    //  if constexpr (std::is_same<S, R>::value) {
    //    if constexpr (std::is_same<T, Q>::value) {
    //        //std::shared_ptr<array::ArrayHandler<S, T>> m_rq = handler;
    //        m_rq = handler;
    //        //std::cout<<"Bla bla"<<std::endl;
    //    } else if constexpr (std::is_same<T, P>::value) {
    //        //std::shared_ptr<array::ArrayHandler<S, T>> m_rp = handler;
    //        m_rp = handler;
    //        //std::cout<<"Bla bla"<<std::endl;
    //    } else if constexpr (std::is_same<T, R>::value) {
    //        //std::shared_ptr<array::ArrayHandler<S, T>> m_rr = handler;
    //        m_rr = handler;
    //        //std::cout<<"Bla bla"<<std::endl;
    //    } else {
    //        error("Container type passed to Builder out of scope of ArrayHandlers container types");
    //    }
    //  } else if constexpr (std::is_same<S, Q>::value) {
    //      if constexpr (std::is_same<T, R>::value) {  // redundant?
    //        //std::shared_ptr<array::ArrayHandler<S, T>> m_qr = handler;
    //        m_qr = handler;
    //        //std::cout<<"Bla bla"<<std::endl;
    //    } else if constexpr (std::is_same<T, P>::value) {
    //        //std::shared_ptr<array::ArrayHandler<S, T>> m_qp = handler;
    //        m_qp = handler;
    //        //std::cout<<"Bla bla"<<std::endl;
    //    } else if constexpr (std::is_same<T, Q>::value) {
    //        //std::shared_ptr<array::ArrayHandler<S, T>> m_qq = handler;
    //        m_qq = handler;
    //    } else {
    //        error("Container type passed to Builder out of scope of ArrayHandlers container types");
    //    }
    //  } else if constexpr (std::is_same<S, P>::value) {
    //      if (std::is_same<T, P>::value) {
    //        //m_pp = std::make_shared<array::ArrayHandler<P, P>>(handler);
    //        m_pp = handler;
    //      std::cout<<"Bla bla"<<std::endl;
    //    } else {
    //        error("With the first container being SparceArray, second can be Sparce Array only");
    //    }
    //  } else {
    //      error("Container type passed to Builder out of scope of ArrayHandlers container types");
    //  }
    //};

    void error(const std::string &message) { throw util::ArrayHandlersError{message}; };
    
  };

  auto& rr() { return *m_rr; } // Do we need the getters?
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
