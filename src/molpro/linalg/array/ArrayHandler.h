#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLER_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLER_H
#include <functional>
#include <list>
#include <memory>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

//! Registers operations identified by their Arguments
//! @tparam Args arguments to the operation. Consider using reference_wrapper if arguments are references.
template <typename... Args> struct OperationRegister {
  using OP = std::tuple<Args...>;
  std::list<std::tuple<Args...>> m_register; //!< ordered register of operations

  //! Push a new operation to the register, collecting all operations marked under Nth argument into consecutive groups.
  //! Those groups are ordered on a first come basis.
  //! @tparam N if N >= 0 order operations by grouping identical Nth arguments together, otherwise no reordering
  //! @tparam ArgEqual functor bool operator()(const Arg&, const Arg&); returning true if two Nth arguments are equal
  template <int N, class ArgEqual> void push(const Args &... args, ArgEqual equal) {
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
template <typename T = int> struct RefEqual {
  bool operator()(const std::reference_wrapper<T> &l, const std::reference_wrapper<T> &r) {
    return std::addressof(l.get()) == std::addressof(r.get());
  }
};

} // namespace util

//! Base- classes and intermediates for ArrayHandler
namespace handler {
template <typename AL, typename AR> class LazyHandleBase;

//! Common base for ArrayHandler with all virtual methods that need to be implemented
template <typename AL, typename AR, bool = !std::is_same<AL, AR>::value> class ArrayHandlerBase {
  friend LazyHandleBase<AL, AR>;

public:
  using value_type_L = typename AL::value_type;
  using value_type_R = typename AR::value_type;
  using value_type = decltype(value_type_L{} * value_type_R{});

protected:
  virtual AL copyLL(const AL &source) = 0;
  virtual AL copyRL(const AR &source) = 0;
  virtual AR copyLR(const AL &source) = 0;
  virtual AR copyRR(const AR &source) = 0;
  virtual AL &scalL(value_type alpha, AL &x) = 0;
  virtual AR &scalR(value_type alpha, AR &x) = 0;
  virtual AL &fillL(value_type alpha, AL &x) = 0;
  virtual AR &fillR(value_type alpha, AR &x) = 0;
  virtual AL &axpyLL(value_type alpha, const AL &x, AL &y) = 0;
  virtual AL &axpyRL(value_type alpha, const AR &x, AL &y) = 0;
  virtual AR &axpyLR(value_type alpha, const AL &x, AR &y) = 0;
  virtual AR &axpyRR(value_type alpha, const AR &x, AR &y) = 0;
  virtual value_type dotLL(const AL &x, const AL &y) = 0;
  virtual value_type dotLR(const AL &x, const AR &y) = 0;
  virtual value_type dotRL(const AR &x, const AL &y) = 0;
  virtual value_type dotRR(const AR &x, const AR &y) = 0;
  virtual void fused_axpyLL(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                            const std::vector<value_type> &alphas,
                            const std::vector<std::reference_wrapper<const AL>> &xx,
                            std::vector<std::reference_wrapper<AL>> &yy) = 0;
  virtual void fused_axpyLR(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                            const std::vector<value_type> &alphas,
                            const std::vector<std::reference_wrapper<const AL>> &xx,
                            std::vector<std::reference_wrapper<AR>> &yy) = 0;
  virtual void fused_axpyRL(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                            const std::vector<value_type> &alphas,
                            const std::vector<std::reference_wrapper<const AR>> &xx,
                            std::vector<std::reference_wrapper<AL>> &yy) = 0;
  virtual void fused_axpyRR(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                            const std::vector<value_type> &alphas,
                            const std::vector<std::reference_wrapper<const AR>> &xx,
                            std::vector<std::reference_wrapper<AR>> &yy) = 0;

  virtual void fused_dotLL(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                           const std::vector<std::reference_wrapper<const AL>> &xx,
                           const std::vector<std::reference_wrapper<const AL>> &yy,
                           std::vector<std::reference_wrapper<value_type>> &out) = 0;
  virtual void fused_dotLR(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                           const std::vector<std::reference_wrapper<const AL>> &xx,
                           const std::vector<std::reference_wrapper<const AR>> &yy,
                           std::vector<std::reference_wrapper<value_type>> &out) = 0;
  virtual void fused_dotRL(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                           const std::vector<std::reference_wrapper<const AR>> &xx,
                           const std::vector<std::reference_wrapper<const AL>> &yy,
                           std::vector<std::reference_wrapper<value_type>> &out) = 0;
  virtual void fused_dotRR(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                           const std::vector<std::reference_wrapper<const AR>> &xx,
                           const std::vector<std::reference_wrapper<const AR>> &yy,
                           std::vector<std::reference_wrapper<value_type>> &out) = 0;
};

//! Interface methods for the ArrayHandler when AL!=AR
template <typename AL, typename AR, bool = !std::is_same<AL, AR>::value>
class ArrayHandlerInterface : public ArrayHandlerBase<AL, AR> {
public:
  using typename ArrayHandlerBase<AL, AR>::value_type_L;
  using typename ArrayHandlerBase<AL, AR>::value_type_R;
  using typename ArrayHandlerBase<AL, AR>::value_type;
  //! Return a copy of an array. This function is only one when AL and AR are the same.
  AL copySame(const AL &source) { return copyLL(source); }
  AL copy(const AR &source) { return copyRL(source); }
  AR copy(const AL &source) { return copyLR(source); }
  AR copySame(const AR &source) { return copyRR(source); }

  //! scale all elements of array x by constant value alpha
  AL &scal(value_type alpha, AL &x) { return scalL(alpha, x); }
  AR &scal(value_type alpha, AR &x) { return scalR(alpha, x); }
  //! replace all elements of array x by constant value alpha
  AL &fill(value_type alpha, AL &x) { return fillL(alpha, x); }
  AR &fill(value_type alpha, AR &x) { return fillR(alpha, x); }
  //! y[:] += alpha * x[:]
  AL &axpy(value_type alpha, const AL &x, AL &y) { return axpyLL(alpha, x, y); }
  AL &axpy(value_type alpha, const AR &x, AL &y) { return axpyRL(alpha, x, y); }
  AR &axpy(value_type alpha, const AL &x, AR &y) { return axpyLR(alpha, x, y); }
  AR &axpy(value_type alpha, const AR &x, AR &y) { return axpyRR(alpha, x, y); }
  //! inner product: result = x.y
  value_type dot(const AL &x, const AL &y) { return dotLL(x, y); }
  value_type dot(const AL &x, const AR &y) { return dotLR(x, y); }
  value_type dot(const AR &x, const AL &y) { return dotRL(x, y); }
  value_type dot(const AR &x, const AR &y) { return dotRR(x, y); }

protected:
  using ArrayHandlerBase<AL, AR>::scalL;
  using ArrayHandlerBase<AL, AR>::scalR;
  using ArrayHandlerBase<AL, AR>::fillL;
  using ArrayHandlerBase<AL, AR>::fillR;
  using ArrayHandlerBase<AL, AR>::copyLL;
  using ArrayHandlerBase<AL, AR>::copyRL;
  using ArrayHandlerBase<AL, AR>::copyLR;
  using ArrayHandlerBase<AL, AR>::copyRR;
  using ArrayHandlerBase<AL, AR>::axpyLL;
  using ArrayHandlerBase<AL, AR>::axpyRL;
  using ArrayHandlerBase<AL, AR>::axpyLR;
  using ArrayHandlerBase<AL, AR>::axpyRR;
  using ArrayHandlerBase<AL, AR>::dotLL;
  using ArrayHandlerBase<AL, AR>::dotLR;
  using ArrayHandlerBase<AL, AR>::dotRL;
  using ArrayHandlerBase<AL, AR>::dotRR;
};

//! Interface methods for the ArrayHandler when AL==AR
template <typename AL, typename AR> class ArrayHandlerInterface<AL, AR, false> : public ArrayHandlerBase<AL, AR> {
public:
  using value_type_L = typename AL::value_type;
  using value_type_R = typename AR::value_type;
  using value_type = decltype(value_type_L{} * value_type_R{});
  //! Return a copy of right array. This function is only one when AL and AR are the same.
  AL copySame(const AL &source) { return copyLL(source); }
  AL copy(const AL &source) { return copySame(source); }

  AL &scal(value_type alpha, AL &x) { return scalL(alpha, x); }
  AL &fill(value_type alpha, AL &x) { return fillL(alpha, x); }
  AL &axpy(value_type alpha, const AL &x, AL &y) { return axpyLL(alpha, x, y); }
  value_type dot(const AL &x, const AL &y) { return dotLL(x, y); }

protected:
  using ArrayHandlerBase<AL, AR>::scalL;
  using ArrayHandlerBase<AL, AR>::fillL;
  using ArrayHandlerBase<AL, AR>::copyLL;
  using ArrayHandlerBase<AL, AR>::axpyLL;
  using ArrayHandlerBase<AL, AR>::dotLL;
};

//! Base class for ArrayHandler::LazyHandle. It contains all the virtual functions and member variables.
template <typename AL, typename AR> class LazyHandleBase {
public:
  using value_type = decltype(typename AL::value_type{} * typename AR::value_type{});
  template <typename T> using ref_wrap = std::reference_wrapper<T>;

  virtual ~LazyHandleBase() = default;

protected:
  //! Types of operations currently registered. Types are strings, because derived classes might add new operations
  std::set<std::string> m_op_types;
  //! register of axpy operations
  util::OperationRegister<value_type, ref_wrap<const AL>, ref_wrap<AL>> m_axpyLL;
  util::OperationRegister<value_type, ref_wrap<const AL>, ref_wrap<AR>> m_axpyLR;
  util::OperationRegister<value_type, ref_wrap<const AR>, ref_wrap<AL>> m_axpyRL;
  util::OperationRegister<value_type, ref_wrap<const AR>, ref_wrap<AR>> m_axpyRR;
  //! register of dot operations
  util::OperationRegister<ref_wrap<const AL>, ref_wrap<const AL>, ref_wrap<value_type>> m_dotLL;
  util::OperationRegister<ref_wrap<const AL>, ref_wrap<const AR>, ref_wrap<value_type>> m_dotLR;
  util::OperationRegister<ref_wrap<const AR>, ref_wrap<const AL>, ref_wrap<value_type>> m_dotRL;
  util::OperationRegister<ref_wrap<const AR>, ref_wrap<const AR>, ref_wrap<value_type>> m_dotRR;

  LazyHandleBase() = default;

  //! Indicates a fatal error occurred
  virtual void error(std::string message) = 0;

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

  virtual void axpyLL(value_type alpha, const AL &x, AL &y) {
    if (register_op_type("axpyLL"))
      m_axpyLL.push(alpha, std::cref(x), std::ref(y));
    else
      error("Failed to register operation type axpyLL with the current state of the LazyHandle");
  }
  virtual void axpyRL(value_type alpha, const AR &x, AL &y) {
    if (register_op_type("axpyRL"))
      m_axpyRL.push(alpha, std::cref(x), std::ref(y));
    else
      error("Failed to register operation type axpyRL with the current state of the LazyHandle");
  }
  virtual void axpyLR(value_type alpha, const AL &x, AR &y) {
    if (register_op_type("axpyLR"))
      m_axpyLR.push(alpha, std::cref(x), std::ref(y));
    else
      error("Failed to register operation type axpyLR with the current state of the LazyHandle");
  }
  virtual void axpyRR(value_type alpha, const AR &x, AR &y) {
    if (register_op_type("axpyRR"))
      m_axpyRR.push(alpha, std::cref(x), std::ref(y));
    else
      error("Failed to register operation type axpyRR with the current state of the LazyHandle");
  }
  virtual void dotLL(const AL &x, const AL &y, value_type &out) {
    if (register_op_type("dotLL"))
      m_dotLL.push(std::cref(x), std::cref(y), std::ref(out));
    else
      error("Failed to register operation type dotLL with the current state of the LazyHandle");
  }
  virtual void dotLR(const AL &x, const AR &y, value_type &out) {
    if (register_op_type("dotLR"))
      m_dotLR.push(std::cref(x), std::cref(y), std::ref(out));
    else
      error("Failed to register operation type dotLR with the current state of the LazyHandle");
  }
  virtual void dotRL(const AR &x, const AL &y, value_type &out) {
    if (register_op_type("dotRL"))
      m_dotRL.push(std::cref(x), std::cref(y), std::ref(out));
    else
      error("Failed to register operation type dotRL with the current state of the LazyHandle");
  }
  virtual void dotRR(const AR &x, const AR &y, value_type &out) {
    if (register_op_type("dotRR"))
      m_dotRR.push(std::cref(x), std::cref(y), std::ref(out));
    else
      error("Failed to register operation type dotRR with the current state of the LazyHandle");
  }

  void clear() {
    m_op_types.clear();
    m_axpyLL.clear();
    m_axpyLR.clear();
    m_axpyRL.clear();
    m_axpyRR.clear();
    m_dotLL.clear();
    m_dotLR.clear();
    m_dotRL.clear();
    m_dotRR.clear();
  }

public:
  virtual void eval(ArrayHandlerBase<AL, AR> &handler) {
    if (!m_axpyLL.empty()) {
      auto reg = util::remove_duplicates<value_type, ref_wrap<const AL>, ref_wrap<AL>, std::equal_to<value_type>,
                                         util::RefEqual<const AL>, util::RefEqual<AL>>(m_axpyLL.m_register, {}, {}, {});
      handler.fused_axpyLL(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
    }
    if (!m_axpyLR.empty()) {
      auto reg = util::remove_duplicates<value_type, ref_wrap<const AL>, ref_wrap<AR>, std::equal_to<value_type>,
                                         util::RefEqual<const AL>, util::RefEqual<AR>>(m_axpyLR.m_register, {}, {}, {});
      handler.fused_axpyLR(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
    }
    if (!m_axpyRL.empty()) {
      auto reg = util::remove_duplicates<value_type, ref_wrap<const AR>, ref_wrap<AL>, std::equal_to<value_type>,
                                         util::RefEqual<const AR>, util::RefEqual<AL>>(m_axpyRL.m_register, {}, {}, {});
      handler.fused_axpyRL(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
    }
    if (!m_axpyRR.empty()) {
      auto reg = util::remove_duplicates<value_type, ref_wrap<const AR>, ref_wrap<AR>, std::equal_to<value_type>,
                                         util::RefEqual<const AR>, util::RefEqual<AR>>(m_axpyRR.m_register, {}, {}, {});
      handler.fused_axpyRR(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
    }
    if (!m_dotLL.empty()) {
      auto reg =
          util::remove_duplicates<ref_wrap<const AL>, ref_wrap<const AL>, ref_wrap<value_type>,
                                  util::RefEqual<const AL>, util::RefEqual<const AL>, util::RefEqual<value_type>>(
              m_dotLL.m_register, {}, {}, {});
      handler.fused_dotLL(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
    }
    if (!m_dotLR.empty()) {
      auto reg =
          util::remove_duplicates<ref_wrap<const AL>, ref_wrap<const AR>, ref_wrap<value_type>,
                                  util::RefEqual<const AL>, util::RefEqual<const AR>, util::RefEqual<value_type>>(
              m_dotLR.m_register, {}, {}, {});
      handler.fused_dotLR(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
    }
    if (!m_dotRL.empty()) {
      auto reg =
          util::remove_duplicates<ref_wrap<const AR>, ref_wrap<const AL>, ref_wrap<value_type>,
                                  util::RefEqual<const AR>, util::RefEqual<const AL>, util::RefEqual<value_type>>(
              m_dotRL.m_register, {}, {}, {});
      handler.fused_dotRL(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
    }
    if (!m_dotRR.empty()) {
      auto reg =
          util::remove_duplicates<ref_wrap<const AR>, ref_wrap<const AR>, ref_wrap<value_type>,
                                  util::RefEqual<const AR>, util::RefEqual<const AR>, util::RefEqual<value_type>>(
              m_dotRR.m_register, {}, {}, {});
      handler.fused_dotRR(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
    }
    clear();
  }
};

//! Interface to ArrayHandler::LazyHandle when AL!=AR.
template <typename AL, typename AR, bool = !std::is_same<AL, AR>::value>
class LazyHandleInterface : public LazyHandleBase<AL, AR> {
protected:
  LazyHandleInterface() = default;

  using LazyHandleBase<AL, AR>::axpyLL;
  using LazyHandleBase<AL, AR>::axpyLR;
  using LazyHandleBase<AL, AR>::axpyRL;
  using LazyHandleBase<AL, AR>::axpyRR;
  using LazyHandleBase<AL, AR>::dotLL;
  using LazyHandleBase<AL, AR>::dotLR;
  using LazyHandleBase<AL, AR>::dotRL;
  using LazyHandleBase<AL, AR>::dotRR;

public:
  using typename LazyHandleBase<AL, AR>::value_type;

  void axpy(value_type alpha, const AL &x, AL &y) { return axpyLL(alpha, x, y); }
  void axpy(value_type alpha, const AL &x, AR &y) { return axpyLR(alpha, x, y); }
  void axpy(value_type alpha, const AR &x, AL &y) { return axpyRL(alpha, x, y); }
  void axpy(value_type alpha, const AR &x, AR &y) { return axpyRR(alpha, x, y); }
  void dot(const AL &x, const AL &y, value_type &out) { return dotLL(x, y, out); }
  void dot(const AL &x, const AR &y, value_type &out) { return dotLR(x, y, out); }
  void dot(const AR &x, const AL &y, value_type &out) { return dotRL(x, y, out); }
  void dot(const AR &x, const AR &y, value_type &out) { return dotRR(x, y, out); }
};

//! Interface to ArrayHandler::LazyHandle when AL==AR.
template <typename AL, typename AR> class LazyHandleInterface<AL, AR, false> : public LazyHandleBase<AL, AR> {
protected:
  LazyHandleInterface() = default;
  using LazyHandleBase<AL, AR>::axpyLL;
  using LazyHandleBase<AL, AR>::dotLL;

public:
  using typename LazyHandleBase<AL, AR>::value_type;

  void axpy(value_type alpha, const AL &x, AL &y) { return axpyLL(alpha, x, y); }
  void dot(const AL &x, const AL &y, value_type &out) { return dotLL(x, y, out); }
};
} // namespace handler

/*!
 * @brief Enhances various operations between pairs of arrays and allows dynamic code injection with uniform interface.
 *
 * The handler interface is symmetric with respect to template parameters, ArrayHandler<AL,AR> should behave the same as
 * ArrayHandler<AR,AL>.
 *
 * Operations between pairs of arrays can often be improved with some hindsight, e.g. lazy evaluation
 * and loop fusing, or might require dynamic information, e.g. file name info for copy constructor of disk array.
 * This class is  not meant primarily to beautify the code, but to allow dynamic code injection without a bloated
 * interface.
 *
 *
 * Use cases
 * ---------
 * Copy constructor
 * @code{.cpp}
 * AL fast;
 * AR slow;
 * ArrayHandler<AL,AR> handler;
 * AR copy_of_fast = handler.copy(fast);
 * @endcode
 *
 * Lazy evaluation of axpy
 * @code{.cpp}
 * std::vector<AL> a1;
 * std::vector<AR> a2;
 * Matrix alpha;
 * // ... Assume initialization was done
 * LazyHandle h = handle.lazy_handle();
 * for (auto i = 0; i < a1.size(); ++i)
 *   for (auto j = 0; j < a2.size(); ++j)
 *     h.axpy(alpha[i,j], a1[i], a2[j]);
 * h.eval(); // or let destructor handle it
 * @endcode
 *
 * Lazy evaluation of dot
 * @code{.cpp}
 * std::vector<AL> a1;
 * std::vector<AR> a2;
 * Matrix overlap;
 * // ... Assume initialization was done
 * LazyHandle h = handle.lazy_handle();
 * for (auto i = 0ul; i < a1.size(); ++i)
 *   for (auto j = 0ul; j < a2.size(); ++j)
 *     h.dot(a1[i], a2[j], overlap[i,j]);
 * h.eval(); // or let the destructor handle it
 * @endcode
 */
template <class AL, class AR = AL> class ArrayHandler : public handler::ArrayHandlerInterface<AL, AR> {
protected:
  ArrayHandler() = default;
  ArrayHandler(const ArrayHandler &) = default;

public:
  using typename handler::ArrayHandlerInterface<AL, AR>::value_type_L;
  using typename handler::ArrayHandlerInterface<AL, AR>::value_type_R;
  using typename handler::ArrayHandlerInterface<AL, AR>::value_type;
  using handler::ArrayHandlerInterface<AL, AR>::copySame;
  using handler::ArrayHandlerInterface<AL, AR>::copy;
  using handler::ArrayHandlerInterface<AL, AR>::scal;
  using handler::ArrayHandlerInterface<AL, AR>::fill;
  using handler::ArrayHandlerInterface<AL, AR>::axpy;
  using handler::ArrayHandlerInterface<AL, AR>::dot;

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
  virtual void error(std::string message) = 0;

  /*!
   * @brief Registers operations for lazy evaluation. Evaluation is triggered by calling eval() or on destruction.
   * @warning The evaluation is performed by the handler, so it must be valid when eval() is called.
   */
  class LazyHandle : public handler::LazyHandleInterface<AL, AR> {
  protected:
    using handler::LazyHandleInterface<AL, AR>::error;

    void error(std::string message) override { m_handler.error(message); };

  public:
    using handler::LazyHandleInterface<AL, AR>::axpy;
    using handler::LazyHandleInterface<AL, AR>::dot;
    LazyHandle(ArrayHandler<AL, AR> &handler) : handler::LazyHandleInterface<AL, AR>(), m_handler{handler} {}
    virtual ~LazyHandle() { LazyHandle::eval(); }

    //! Calls handler to evaluate the registered operations
    virtual void eval() {
      if (m_invalid)
        return;
      handler::LazyHandleInterface<AL, AR>::eval(m_handler);
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

    template <typename... Args> void axpy(Args &&... args) {
      m_lazy_handle->axpy(std::forward<Args>(args)...);
      if (m_off)
        eval();
    }
    template <typename... Args> void dot(Args &&... args) {
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

public:
  //! Returns a lazy handle
  virtual ProxyHandle lazy_handle() = 0;
};

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLER_H
