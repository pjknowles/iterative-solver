#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLER_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLER_H
#include <functional>
#include <list>
#include <memory>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

//! Places a new operation in the register, so that high priority arrays that are equal are grouped together
template <typename Op, int N> struct RegisterOperation {
  void operator()(std::list<Op> &reg, Op &&new_op) {
    auto ref = std::get<N>(new_op);
    auto rend = std::find_if(reg.rbegin(), reg.rend(), [&ref](const auto &el) {
      auto xx = std::get<N>(el);
      return std::addressof(ref.get()) == std::addressof(xx.get());
    });
    auto end_of_group = reg.end();
    if (rend != reg.rend())
      end_of_group = rend.base();
    reg.insert(end_of_group, std::forward<Op>(new_op));
  }
};

template <typename Op> struct RegisterOperation<Op, -1> {
  void operator()(std::list<Op> &reg, Op &&new_op) { reg.insert(reg.end(), std::forward<Op>(new_op)); }
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
template <typename X, typename Y, typename S>
std::tuple<std::vector<std::pair<size_t, size_t>>, std::vector<X>, std::vector<Y>, std::vector<S>>
remove_duplicates(const std::list<std::tuple<X, Y, S>> &reg) {
  auto n_op = reg.size();
  std::vector<S> scalar;
  scalar.reserve(n_op);
  std::transform(cbegin(reg), cend(reg), std::back_inserter(scalar), [](const auto &el) { return std::get<2>(el); });
  std::vector<std::pair<size_t, size_t>> op_register;
  op_register.reserve(n_op);
  std::vector<X> xx;
  std::vector<Y> yy;
  for (const auto &op : reg) {
    auto x = std::get<0>(op);
    auto y = std::get<1>(op);
    auto it_x = std::find_if(cbegin(xx), cend(xx),
                             [&x](const auto &el) { return std::addressof(x.get()) == std::addressof(el.get()); });
    auto it_y = std::find_if(cbegin(yy), cend(yy),
                             [&y](const auto &el) { return std::addressof(y.get()) == std::addressof(el.get()); });
    auto ix = distance(cbegin(xx), it_x);
    auto iy = distance(cbegin(yy), it_y);
    if (it_x == cend(xx))
      xx.push_back(x);
    if (it_y == cend(yy))
      yy.push_back(y);
    op_register.emplace_back(ix, iy);
  }
  return {op_register, xx, yy, scalar};
}

} // namespace util

/*!
 * @brief Enhances various operations between pairs of arrays and allows dynamic code injection with uniform interface.
 *
 * Operations between pairs of arrays can often be improved with some hindsight, e.g. lazy evaluation
 * and loop fusing, or might require dynamic information, e.g. file name info for copy constructor of disk array.
 *
 * This class is designed specifically to deal with operations that require pairs of arrays. It is not meant primarily
 * to beautify the code, but to allow dynamic code injection without a bloated interface.
 *
 * Concepts
 * --------
 *   * Left array (AL) - array to the left of a binary operation
 *   * Right array (AR) - array to the right of a binary operation
 *
 * For example:
 * @code{.cpp}
 * //  x[:] += a * y[:]
 * void axpy(LeftArray x, value_type a, RightArray y);
 * // copy constructor
 * LeftArray(const RightArray&);
 * @endcode
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
 * for (auto i = 0ul; i < a1.size(); ++i)
 *   for (auto j = 0ul; j < a2.size(); ++j)
 *     h.axpy(a1[i], alpha[i,j], a2[j]);
 * // h.dot(a1[0], a2[0], alpha[0]); // <- should throw an error, cannot mix operations. Just start a new handle!
 * h.eval(); // or let destructor handle it
 * @endcode
 *
 * Lazy evaluation of dot
 * @code{.cpp}
 * std::vector<AL> a1;
 * std::vector<AR> a2;
 * Matrix result;
 * // ... Assume initialization was done
 * LazyHandle h = handle.lazy_handle();
 * for (auto i = 0ul; i < a1.size(); ++i)
 *   for (auto j = 0ul; j < a2.size(); ++j)
 *     h.dot(a1[i], a2[j], result[i,j]);
 * h.eval(); // or let the destructor handle it
 * @endcode
 *
 * Design
 * ------
 * The numerical work should be done by the handler.
 * Lazy evaluation operations can be registered with a LazyHandle.
 * It must be passed as a pointer to allow inheritance. Returning it as reference will trigger copy constructor when
 * used with auto.
 *
 *
 * @tparam AL left array type
 * @tparam AR right array type
 */
template <class AL, class AR = AL> class ArrayHandler {

public:
  using value_type = typename AL::value_type;
  ArrayHandler(const ArrayHandler &) = delete;

protected:
  ArrayHandler() = default;

public:
  virtual ~ArrayHandler() {
    std::for_each(m_lazy_handles.begin(), m_lazy_handles.end(), [](auto &el) {
      if (auto handle = el.lock())
        handle->invalidate();
    });
  }

  //! Return a copy of right array. This function is only one when AL and AR are the same.
  template <bool flag = std::is_same<AL, AR>::value, typename = typename std::enable_if_t<flag>>
  AR copy(const AR &source) {
    return copySame(source);
  }
  //! Return a copy of right array as a left array. This function is only one when AL and AR are different.
  template <bool flag = std::is_same<AL, AR>::value, typename std::enable_if_t<!flag, nullptr_t> = nullptr>
  AL copy(const AR &source) {
    return copyR(source);
  }
  //! Return a copy of left array as a right array. This function is only one when AL and AR are different.
  template <bool flag = !std::is_same<AL, AR>::value, typename std::enable_if_t<flag, nullptr_t> = nullptr>
  AR copy(const AL &source) {
    return copyL(source);
  }

  virtual AL &axpy(AL &x, value_type alpha, const AR &y) = 0;
  virtual value_type dot(const AL &x, const AR &y) = 0;

protected:
  virtual AR copySame(const AR &source) = 0;
  virtual AL copyR(const AR &source) = 0;
  virtual AR copyL(const AL &source) = 0;
  /*!
   * @brief
   * @param reg register of unique operations. For each element first is the index to x in xx, second is index to y in
   * yy
   * @param xx references to unique x arrays
   * @param alphas full list of alphas
   * @param yy references to unique y arrays
   */
  virtual void fused_axpy(const std::vector<std::pair<size_t, size_t>> &reg,
                          std::vector<std::reference_wrapper<AL>> &xx,
                          std::vector<std::reference_wrapper<const AR>> &yy, std::vector<value_type> alphas) = 0;

  virtual void fused_dot(const std::vector<std::pair<size_t, size_t>> &reg,
                         std::vector<std::reference_wrapper<const AL>> &xx,
                         std::vector<std::reference_wrapper<const AR>> &yy,
                         std::vector<std::reference_wrapper<value_type>> &out) = 0;

  /*!
   * @brief Throws an error
   * @param message error message
   */
  virtual void error(std::string message) = 0;

  //! Assigns highest priority for accessing local buffer to either left or right array, or non to keep the original
  //! order of registered operations
  enum class Priority { left, right, none };

  //! Registers operations for lazy evaluation. Evaluation is triggered by calling eval() or on destruction.
  //! @warning The evaluation is performed by the handler, so it must be valid when eval() is called.
  class LazyHandle {
  protected:
    using OPaxpy = std::tuple<std::reference_wrapper<AL>, std::reference_wrapper<const AR>, value_type>;
    using OPdot = std::tuple<std::reference_wrapper<const AL>, std::reference_wrapper<const AR>,
                             std::reference_wrapper<value_type>>;

  public:
    LazyHandle(ArrayHandler<AL, AR> &handler, Priority priority) : m_handler{handler}, m_priority{priority} {}
    virtual ~LazyHandle() { LazyHandle::eval(); }

    //! Registers an axpy() operation
    virtual void axpy(AL &x, value_type alpha, const AR &y) {
      if (m_invalid)
        return;
      if (!m_dot_register.empty())
        m_handler.error("Trying to register axpy when dot is already assigned. LazyHandle can only register "
                        "operations of the same type.");
      if (m_priority == Priority::left)
        util::RegisterOperation<OPaxpy, 0>()(m_axpy_register, {x, y, alpha});
      else if (m_priority == Priority::right)
        util::RegisterOperation<OPaxpy, 1>()(m_axpy_register, {x, y, alpha});
      else
        util::RegisterOperation<OPaxpy, -1>()(m_axpy_register, {x, y, alpha});
    }

    //! Registers a dot() operation
    virtual void dot(const AL &x, const AR &y, value_type &out) {
      if (m_invalid)
        return;
      if (!m_axpy_register.empty())
        m_handler.error("Trying to register dot when axpy is already assigned. LazyHandle can only register "
                        "operations of the same type.");
      if (m_priority == Priority::left)
        util::RegisterOperation<OPdot, 0>()(m_dot_register, {x, y, std::ref(out)});
      else if (m_priority == Priority::right)
        util::RegisterOperation<OPdot, 1>()(m_dot_register, {x, y, std::ref(out)});
      else
        util::RegisterOperation<OPdot, -1>()(m_dot_register, {x, y, std::ref(out)});
    }

    //! Calls handler to evaluate the registered operations
    virtual void eval() {
      if (m_invalid)
        return;
      if (!m_axpy_register.empty()) {
        auto reg = util::remove_duplicates(m_axpy_register);
        m_handler.fused_axpy(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
      }
      if (!m_dot_register.empty()) {
        auto reg = util::remove_duplicates(m_dot_register);
        m_handler.fused_dot(std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg));
      }
    }

    //! Flag the handler as invalid so that no new operations are registered operations eval() does nothing
    void invalidate() { m_invalid = true; }
    //! Returns true if the handle is marked as invalid. LazyHandle becomes invalid when the overlying ArrayHandler is
    //! destroyed, or invalidate() is called
    bool invalid() { return m_invalid; }

  protected:
    ArrayHandler<AL, AR> &m_handler;      //!< all operations are still done through the handler
    bool m_invalid = false;               //!< flags if the handler has been destroyed and LazyHandle is now invalid
    Priority m_priority = Priority::none; //!< how to prioritise evaluation when fusing a loop
    std::list<OPaxpy> m_axpy_register;    //!< register of axpy operations
    std::list<OPdot> m_dot_register;      //!< register of dot operations
  };

  //! A convenience wrapper around a pointer to the LazyHandle
  class ProxyHandle {
  public:
    ProxyHandle(std::shared_ptr<LazyHandle> handle) : m_lazy_handle{std::move(handle)} {}

    template <typename... Args> void axpy(Args &&... args) { m_lazy_handle->axpy(std::forward<Args>(args)...); }
    template <typename... Args> void dot(Args &&... args) { m_lazy_handle->dot(std::forward<Args>(args)...); }
    void eval() { m_lazy_handle->eval(); }
    void invalidate() { m_lazy_handle->invalidate(); }
    bool invalid() { return m_lazy_handle->invalid(); }

  protected:
    std::shared_ptr<LazyHandle> m_lazy_handle;
  };

  std::vector<std::weak_ptr<LazyHandle>>
      m_lazy_handles; //!< keeps track of all created lazy handles so they can be invalidated on destruction.

public:
  //! Returns a lazy handle
  virtual ProxyHandle lazy_handle() {
    auto handle = std::make_shared<LazyHandle>(*this, Priority::none);
    auto empty_handle =
        std::find_if(m_lazy_handles.begin(), m_lazy_handles.end(), [](const auto &el) { return el.expired(); });
    if (empty_handle == m_lazy_handles.end())
      m_lazy_handles.push_back(handle);
    else
      *empty_handle = handle;
    return handle;
  };
};

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLER_H
