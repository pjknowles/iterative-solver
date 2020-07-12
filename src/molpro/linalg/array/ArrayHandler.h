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
template <typename Op, int N> static void register_op(std::list<Op> &reg, Op &&new_op) {
  auto end_of_group = reg.end();
  if (N >= 0) {
    auto ref = std::get<N>(new_op);
    end_of_group = std::find_if(reg.rbegin(), reg.rend(), [&ref](const auto &el) {
      auto xx = std::get<N>(el);
      return std::addressof(ref.get()) == std::addressof(xx.get());
    });
  }
  if (end_of_group == reg.rbegin())
    end_of_group = reg.end();
  reg.insert(end_of_group, std::forward<Op>(new_op));
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
 * //OR could use a proxy
 * copy_of_fast = handler(fast);
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
 *     h.dot(a1[i], alpha[i,j], a2[j]);
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
template <class AL, class AR> class ArrayHandler {
public:
  ArrayHandler(const ArrayHandler &) = delete;
  virtual ~ArrayHandler() {
    std::transform(m_lazy_handles.begin(), m_lazy_handles.end(), [](auto el) {
      if (auto handle = el.lock())
        handle->invalidate();
    });
  }

  virtual AL copy(const AR &source) = 0;
  virtual AR copy(const AL &source) = 0;

  virtual AL &axpy(typename AR::value_type alpha, const AR &y) = 0;
  virtual AL &dot(typename AR::value_type alpha, const AR &y) = 0;

protected:
  /*!
   * @brief
   * @param reg register of unique operations. For each element first is the index to x in xx, second is index to y in
   * yy
   * @param xx references to unique x arrays
   * @param alphas full list of alphas
   * @param yy references to unique y arrays
   */
  virtual void fused_axpy(const std::vector<std::tuple<size_t, size_t>> &reg,
                          std::vector<std::reference_wrapper<AL>> &xx, std::vector<typename AR::value_type> alphas,
                          std::vector<std::reference_wrapper<const AR>> &yy) = 0;

  virtual void fused_dot(const std::vector<std::tuple<size_t, size_t>> &reg,
                         std::vector<std::reference_wrapper<AL>> &xx, std::vector<std::reference_wrapper<const AR>> &yy,
                         std::vector<std::reference_wrapper<typename AR::value_type>> &out) = 0;

  //! Assigns highest priority for accessing local buffer to either left or right array, or non to keep the original
  //! order of registered operations
  enum class Priority { left, right, none };

  //! Registers operations for lazy evaluation. Evaluation is triggered by calling eval() or on destruction.
  //! @warning The evaluation is performed by the handler, so it must be valid when eval() is called.
  class LazyHandle {
  public:
    LazyHandle(const ArrayHandler<AL, AR> &handler, Priority priority) : m_handler{handler}, m_priority{priority} {}
    virtual ~LazyHandle() { LazyHandle::eval(); }

    //! Registers an axpy() operation
    virtual void axpy(AL &x, typename AR::value_type alpha, const AR &y) {
      if (m_invalid)
        return;
      if (!m_dot_register.empty())
        m_handler.raise_error("Trying to register axpy when dot is already assigned. LazyHandle can only register "
                              "operations of the same type.");
      if (m_priority == Priority::left)
        util::register_op<decltype(m_axpy_register.front()), 0>(m_axpy_register, {x, y, alpha});
      else if (m_priority == Priority::right)
        util::register_op<decltype(m_axpy_register.front()), 1>(m_axpy_register, {x, y, alpha});
      else
        util::register_op<decltype(m_axpy_register.front()), -1>(m_axpy_register, {x, y, alpha});
    }

    //! Registers a dot() operation
    virtual void dot(const AL &x, const AR &y, typename AR::value_type &out) {
      if (m_invalid)
        return;
      if (!m_axpy_register.empty())
        m_handler.raise_error("Trying to register dot when axpy is already assigned. LazyHandle can only register "
                              "operations of the same type.");
      if (m_priority == Priority::left)
        util::register_op<decltype(m_dot_register.front()), 0>(m_dot_register, {x, y, std::ref(out)});
      else if (m_priority == Priority::right)
        util::register_op<decltype(m_dot_register.front()), 1>(m_dot_register, {x, y, std::ref(out)});
      else
        util::register_op<decltype(m_dot_register.front()), -1>(m_dot_register, {x, y, std::ref(out)});
    }

    //! Sorts the registered operations to prioritize access to the left or right array, or keep in the original order,
    //! merges operations if possible (e.g. same x, y, different alpha) and calls handler to evaluate
    virtual void eval() {
      if (m_invalid)
        return;
      if (!m_axpy_register.empty()) {
        // Find duplicates and store unique elements in a separate vector
        std::vector<std::reference_wrapper<typename AR::value_type>> scalar;
        scalar.reserve(m_axpy_register.size());
        std::copy(m_axpy_register.cbegin(), m_axpy_register.cend(), std::back_inserter(scalar),
                  [](const auto &el) { return std::get<2>(el); });
        std::vector<std::pair<size_t, size_t>> op_register; // pair of x and y array indices for each operation
        op_register.resize(m_axpy_register.size());
        std::vector<std::reference_wrapper<AL>> xx;
        std::vector<std::reference_wrapper<const AR>> yy;
        for (size_t i = 0; i < m_axpy_register.size(); ++i) {
          auto x = std::get<0>(m_axpy_register[i]);
          auto y = std::get<1>(m_axpy_register[i]);
          auto it_x = std::find_if(cbegin(xx), cend(xx), [&x](const auto &el) {
            return std::addressof(x.get()) == std::addressof(el.get());
          });
          auto it_y = std::find_if(cbegin(yy), cend(yy), [&y](const auto &el) {
            return std::addressof(y.get()) == std::addressof(el.get());
          });
          auto ix = distance(cbegin(xx), it_x);
          auto iy = distance(cbegin(yy), it_y);
          if (it_x == cend(xx))
            xx.push_back(x);
          if (it_y == cend(yy))
            yy.push_back(y);
          op_register.emplace_back(ix, iy);
        }
        m_handler.fused_axpy(op_register, xx, yy, scalar);
      }
      if (!m_dot_register.empty())
        m_handler.fused_dot();
    }

    //! Invalidate the handler so that registered operations do not get evaluated
    void invalidate() { m_invalid = true; }
    //! Returns true if the handle is invalid. LazyHandle becomes invalid when the overlying ArrayHandler is destroyed,
    //! or invalidate() is called
    bool invalid() { return m_invalid; }

  protected:
    ArrayHandler<AL, AR> &m_handler;      //! all operations are still done through the handler
    bool m_invalid = false;               //! flags if the handler has been destroyed and LazyHandle is now invalid
    Priority m_priority = Priority::none; //! how to prioritise evaluation when fusing a loop
    std::list<std::tuple<std::reference_wrapper<AL>, std::reference_wrapper<const AR>, typename AR::value_type>>
        m_axpy_register; //! register of axpy operations
    std::list<std::tuple<std::reference_wrapper<const AL>, std::reference_wrapper<const AR>,
                         std::reference_wrapper<typename AR::value_type>>>
        m_dot_register; //! register of dot operations
  };

  class ProxyHandle {
  protected:
    ProxyHandle(std::shared_ptr<LazyHandle> handle) : m_lazy_handle{std::move(handle)} {}

  public:
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
