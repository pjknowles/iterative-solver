#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLER_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLER_H
#include <algorithm>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <numeric>
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif

#include <molpro/linalg/array/type_traits.h>
#include <molpro/linalg/itsolv/wrap_util.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>

using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::subspace::Matrix;

namespace molpro::linalg::array {
namespace util {

struct ArrayHandlerError : public std::logic_error {
  using std::logic_error::logic_error;
};

//! Registers operations identified by their Arguments
//! @tparam Args arguments to the operation. Consider using reference_wrapper if arguments are references.
template <typename... Args>
struct OperationRegister {
  using OP = std::tuple<Args...>;
  std::list<std::tuple<Args...>> m_register; //!< ordered register of operations

  //! Push a new operation to the register, collecting all operations marked under Nth argument into consecutive groups.
  //! Those groups are ordered on a first come basis.
  //! @tparam N if N >= 0 order operations by grouping identical Nth arguments together, otherwise no reordering
  //! @tparam ArgEqual functor bool operator()(const Arg&, const Arg&); returning true if two Nth arguments are equal
  template <int N, class ArgEqual>
  void push(const Args &... args, ArgEqual equal) {
    auto &&new_op = OP{args...};
    auto &ref = std::get<N>(new_op);
    auto rend = std::find_if(m_register.rbegin(), m_register.rend(), [&ref, &equal](const auto &el) {
      auto xx = std::get<N>(el);
      return equal(ref, xx);
    });
    auto end_of_group = m_register.end();
    if (rend != m_register.rend())
      end_of_group = rend.base();
    m_register.insert(end_of_group, std::forward<OP>(new_op));
  }

  //! Register each operation as it comes with no reordering.
  void push(const Args &... args) { m_register.push_back({args...}); }

  bool empty() { return m_register.empty(); }
  void clear() { m_register.clear(); }
};

/*!
 * @brief Find duplicates references to x and y arrays and store unique elements in a separate vector.
 *
 * Vector of scalars is returned complete without truncation.
 * The register becomes a vector of pairs, first element of a pair containing index of x array in the vector of unique x
 * arrays, the second is for y arrays.
 *
 * @return Returns the new register in terms of x and y indices, vectors of unique x and y arrays, and vector of all
 * scalars.
 */
//!
template <typename X, typename Y, typename Z, class EqualX, class EqualY, class EqualZ>
std::tuple<std::vector<std::tuple<size_t, size_t, size_t>>, std::vector<X>, std::vector<Y>, std::vector<Z>>
remove_duplicates(const std::list<std::tuple<X, Y, Z>> &reg, EqualX equal_x, EqualY equal_y, EqualZ equal_z) {
  auto n_op = reg.size();
  std::vector<std::tuple<size_t, size_t, size_t>> op_register;
  op_register.reserve(n_op);
  std::vector<X> xx;
  std::vector<Y> yy;
  std::vector<Z> zz;
  for (const auto &op : reg) {
    auto x = std::get<0>(op);
    auto y = std::get<1>(op);
    auto z = std::get<2>(op);
    auto it_x = std::find_if(cbegin(xx), cend(xx), [&x, &equal_x](const auto &el) { return equal_x(x, el); });
    auto it_y = std::find_if(cbegin(yy), cend(yy), [&y, &equal_y](const auto &el) { return equal_y(y, el); });
    auto it_z = std::find_if(cbegin(zz), cend(zz), [&z, &equal_z](const auto &el) { return equal_z(z, el); });
    auto ix = distance(cbegin(xx), it_x);
    auto iy = distance(cbegin(yy), it_y);
    auto iz = distance(cbegin(zz), it_z);
    if (it_x == cend(xx))
      xx.push_back(x);
    if (it_y == cend(yy))
      yy.push_back(y);
    if (it_z == cend(zz))
      zz.push_back(z);
    op_register.emplace_back(ix, iy, iz);
  }
  return {op_register, xx, yy, zz};
}

//! When called returns true if addresses of two references are the same
template <typename T = int>
struct RefEqual {
  bool operator()(const std::reference_wrapper<T> &l, const std::reference_wrapper<T> &r) {
    return std::addressof(l.get()) == std::addressof(r.get());
  }
};
} // namespace util

/*!
 * @brief Enhances various operations between pairs of arrays and allows dynamic code injection with uniform interface.
 *
 * The handler is directional, there is a left and a right array.
 *
 * @tparam AL
 * @tparam AR
 *
 * Operations between pairs of arrays can often be improved with some hindsight, e.g. lazy evaluation
 * and loop fusing, or might require dynamic information, e.g. file name info for copy constructor of disk array.
 * This class is  not meant to beautify the code, but to allow dynamic code injection without a bloated interface.
 *
 *
 *
 * Use cases
 * ---------
 * Copy constructor
 * @code{.cpp}
 * AL slow;
 * AR fast;
 * ArrayHandler<AL, AR> handler;
 * AL copy_of_fast = handler.copy(fast);
 * @endcode
 *
 * Lazy evaluation of axpy
 * @code{.cpp}
 * std::vector<AL> a1;
 * std::vector<AR> a2;
 * Matrix alpha;
 * // ... Assume initialization was done
 * auto h = handle.lazy_handle();
 * for (auto i = 0; i < a1.size(); ++i)
 *   for (auto j = 0; j < a2.size(); ++j)
 *     h.axpy(alpha[i,j], a2[i], a1[j]);
 * h.eval(); // or let destructor handle it
 * @endcode
 *
 * Lazy evaluation of dot
 * @code{.cpp}
 * std::vector<AL> a1;
 * std::vector<AR> a2;
 * Matrix overlap;
 * // ... Assume initialization was done
 * auto h = handle.lazy_handle();
 * for (auto i = 0; i < a1.size(); ++i)
 *   for (auto j = 0; j < a2.size(); ++j)
 *     h.dot(a1[i], a2[j], overlap[i,j]);
 * h.eval(); // or let the destructor handle it
 * @endcode
 */
template <class AL, class AR = AL>
class ArrayHandler {
protected:
  //ArrayHandler() = default;
  ArrayHandler() : m_counter(std::make_shared<Counter>()) {};
  ArrayHandler(const ArrayHandler &) = default;
  
  struct Counter {
    int scal = 0;
    int dot = 0;
    int axpy = 0;
    int copy = 0;
    int gemm_inner = 0;
    int gemm_outer = 0;
  };
  
  std::shared_ptr<Counter> m_counter;
public:
  using value_type_L = typename array::mapped_or_value_type_t<AL>;
  using value_type_R = typename array::mapped_or_value_type_t<AR>;
  using value_type = decltype(value_type_L{} * value_type_R{});
  using value_type_abs = decltype(check_abs<value_type>());

  virtual AL copy(const AR &source) = 0;
  //! Copy content of y into x
  virtual void copy(AL &x, const AR &y) = 0;
  virtual void scal(value_type alpha, AL &x) = 0;
  virtual void fill(value_type alpha, AL &x) = 0;
  virtual void axpy(value_type alpha, const AR &x, AL &y) = 0;
  virtual value_type dot(const AL &x, const AR &y) = 0;
  
  /*!
   * Perform axpy() on multiple pairs of containers in an efficient manner
   */
  virtual void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) = 0;
  
  /*!
   * Perform dot() on multiple pairs of containers in an efficient manner
   */
  virtual Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) = 0;
  
  /*!
   * @brief Select n indices with largest by absolute value contributions to the dot product
   *
   * This is necessary for perturbation theory.
   *
   * @param n number of indices to select
   * @param x left array
   * @param y right array
   * @return map of indices and corresponding x,y product
   */
  virtual std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) = 0;
  
  Counter& counter() {return *m_counter;}

  //! Destroys ArrayHandler instance and invalidates any LazyHandler it created. Invalidated handler will not evaluate.
  virtual ~ArrayHandler() {
    std::for_each(m_lazy_handles.begin(), m_lazy_handles.end(), [](auto &el) {
      if (auto handle = el.lock())
        handle->invalidate();
    });
  }

protected:
  /*!
   * @brief Throws an error
   * @param message error message
   */
  virtual void error(const std::string &message) { throw util::ArrayHandlerError{message}; };

  //! Default implementation of fused_axpy without any simplification
  virtual void fused_axpy(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                          const std::vector<value_type> &alphas,
                          const std::vector<std::reference_wrapper<const AR>> &xx,
                          std::vector<std::reference_wrapper<AL>> &yy) {
    for (const auto &i : reg) {
      size_t ai, xi, yi;
      std::tie(ai, xi, yi) = i;
      axpy(alphas[ai], xx[xi].get(), yy[yi].get());
    }
  }

  //! Default implementation of fused_dot without any simplification
  virtual void fused_dot(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                         const std::vector<std::reference_wrapper<const AL>> &xx,
                         const std::vector<std::reference_wrapper<const AR>> &yy,
                         std::vector<std::reference_wrapper<value_type>> &out) {
    for (const auto &i : reg) {
      size_t xi, yi, zi;
      std::tie(xi, yi, zi) = i;
      out[zi].get() = dot(xx[xi].get(), yy[yi].get());
    }
  }
  

  /*!
   * @brief Registers operations for lazy evaluation. Evaluation is triggered by calling eval() or on destruction.
   * @warning The evaluation is performed by the handler, so it must be valid when eval() is called.
   */
  class LazyHandle {
  public:
    using value_type = ArrayHandler<AL, AR>::value_type;
    template <typename T>
    using ref_wrap = std::reference_wrapper<T>;

  protected:
    //! Types of operations currently registered. Types are strings, because derived classes might add new operations
    std::set<std::string> m_op_types;
    //! register of axpy operations
    util::OperationRegister<value_type, ref_wrap<const AR>, ref_wrap<AL>> m_axpy;
    //! register of dot operations
    util::OperationRegister<ref_wrap<const AL>, ref_wrap<const AR>, ref_wrap<value_type>> m_dot;

    void error(std::string message) { m_handler.error(message); };

    /*!
     * @brief Register an operation type
     *
     * Currently, only one type of operation is suported at any time.
     *
     * @param type type of operation
     * @returns true if that operation type is currently allowed, false otherwise
     */
    virtual bool register_op_type(const std::string &type) {
      if (m_op_types.count(type) == 0 && !m_op_types.empty())
        return false;
      m_op_types.insert(type);
      return true;
    }

    //! Clear the registry
    void clear() {
      m_op_types.clear();
      m_axpy.clear();
      m_dot.clear();
    }

  public:
    explicit LazyHandle(ArrayHandler<AL, AR> &handler) : m_handler{handler} {}
    virtual ~LazyHandle() { LazyHandle::eval(); }

    virtual void axpy(value_type alpha, const AR &x, AL &y) {
      if (register_op_type("axpy"))
        m_axpy.push(alpha, std::cref(x), std::ref(y));
      else
        error("Failed to register operation type axpy with the current state of the LazyHandle");
    }
    virtual void dot(const AL &x, const AR &y, value_type &out) {
      if (register_op_type("dotLR"))
        m_dot.push(std::cref(x), std::cref(y), std::ref(out));
      else
        error("Failed to register operation type dot with the current state of the LazyHandle");
    }

    //! Calls handler to evaluate the registered operations
    virtual void eval() {
      if (m_invalid)
        return;
      if (!m_axpy.empty()) {
        auto reg = util::remove_duplicates<value_type, ref_wrap<const AR>, ref_wrap<AL>, std::equal_to<value_type>,
                                           util::RefEqual<const AR>, util::RefEqual<AL>>(m_axpy.m_register, {}, {}, {});
        m_handler.fused_axpy(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
      }
      if (!m_dot.empty()) {
        auto reg =
            util::remove_duplicates<ref_wrap<const AL>, ref_wrap<const AR>, ref_wrap<value_type>,
                                    util::RefEqual<const AL>, util::RefEqual<const AR>, util::RefEqual<value_type>>(
                m_dot.m_register, {}, {}, {});
        m_handler.fused_dot(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
      }
      clear();
    }

    //! Flag the handler as invalid so that no new operations are registered operations eval() does nothing
    void invalidate() { m_invalid = true; }
    //! Returns true if the handle is marked as invalid. LazyHandle becomes invalid when the overlying ArrayHandler is
    //! destroyed, or invalidate() is called
    bool invalid() { return m_invalid; }

  protected:
    ArrayHandler<AL, AR> &m_handler; //!< all operations are still done through the handler
    bool m_invalid = false;          //!< flags if the handler has been destroyed and LazyHandle is now invalid
  };

  //! A convenience wrapper around a pointer to the LazyHandle
  class ProxyHandle {
  public:
    ProxyHandle(std::shared_ptr<LazyHandle> handle) : m_lazy_handle{std::move(handle)} {}

    template <typename... Args>
    void axpy(Args &&... args) {
      m_lazy_handle->axpy(std::forward<Args>(args)...);
      if (m_off)
        eval();
    }
    template <typename... Args>
    void dot(Args &&... args) {
      m_lazy_handle->dot(std::forward<Args>(args)...);
      if (m_off)
        eval();
    }
    void eval() { m_lazy_handle->eval(); }
    void invalidate() { m_lazy_handle->invalidate(); }
    bool invalid() { return m_lazy_handle->invalid(); }

    //! Turn off lazy evaluation. Next operation will evaluate without delay.
    void off() { m_off = true; };
    //! Turn on lazy evaluation
    void on() { m_off = false; };
    //! Returns true if lazy evaluation is off
    bool is_off() { return m_off; }

  protected:
    std::shared_ptr<LazyHandle> m_lazy_handle;
    bool m_off = false; //!< whether lazy evaluation is on or off
  };

  std::vector<std::weak_ptr<LazyHandle>> m_lazy_handles; //!< keeps track of all created lazy handles

  //! Save  weak ptr to a lazy handle
  void save_handle(const std::shared_ptr<LazyHandle> &handle) {
    auto empty_handle =
        std::find_if(m_lazy_handles.begin(), m_lazy_handles.end(), [](const auto &el) { return el.expired(); });
    if (empty_handle == m_lazy_handles.end())
      m_lazy_handles.push_back(handle);
    else
      *empty_handle = handle;
  }

  ProxyHandle lazy_handle(ArrayHandler<AL, AR> &handler) {
    auto handle = std::make_shared<typename ArrayHandler<AL, AR>::LazyHandle>(handler);
    save_handle(handle);
    return handle;
  };

public:
  //! Returns a lazy handle. Most implementations simply need to call the overload: return lazy_handle(*this);.
  virtual ProxyHandle lazy_handle() = 0;
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLER_H
