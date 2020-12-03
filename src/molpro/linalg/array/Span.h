#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_SPAN_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_SPAN_H
#include <stddef.h>
#include <type_traits>
#include <utility>

namespace molpro::linalg::array {

#if __cplusplus >= 202002L
#include <span>

//! For those who moved on to c++20
using Span = std::span;

#endif

#if __cplusplus < 202002L
inline
#endif
    namespace span {

/*!
 * @brief Non-owning container taking a pointer to the data buffer and its size and exposing routines for iteration
 *
 * This is a poor mimic of std::span. If c++20 standard is used than std::span is aliased to Span.
 */
template <typename T = double>
class Span {
public:
  using element_type = T;
  using value_type = std::remove_cv_t<T>;
  using reference = T&;
  using const_reference = T&;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using iterator = T*;
  using const_iterator = T const*;

  Span() = default;
  ~Span() = default;
  Span(T* data, size_type size) : m_buffer{data}, m_size{size} {}
  Span(const Span& source) = default;
  Span(Span<T>&& source) noexcept : m_buffer{source.m_buffer}, m_size{source.m_size} {
    Span<T> t{};
    swap(source, t);
  }
  Span& operator=(const Span& source) = default;
  Span& operator=(Span&& source) noexcept {
    Span<T> t{std::move(source)};
    swap(*this, t);
    return *this;
  }
  reference operator[](size_type i) { return *(m_buffer + i); }
  const_reference operator[](size_type i) const { return *(m_buffer + i); }

  //! Swap content of two Spans
  friend void swap(Span<T>& x, Span<T>& y) {
    using std::swap;
    swap(x.m_buffer, y.m_buffer);
    swap(x.m_size, y.m_size);
  }

  iterator data() { return m_buffer; }
  const_iterator data() const { return m_buffer; }

  iterator begin() { return m_buffer; }
  const_iterator begin() const { return m_buffer; }
  const_iterator cbegin() const { return m_buffer; }

  iterator end() { return m_buffer + m_size; }
  const_iterator end() const { return m_buffer + m_size; }
  const_iterator cend() const { return m_buffer + m_size; }

  size_type size() const { return m_size; }

  bool empty() const { return m_size == 0; }

protected:
  iterator m_buffer = nullptr;
  size_type m_size = 0;
};

template <typename T>
auto begin(Span<T>& x) {
  return x.begin();
}

template <typename T>
auto begin(const Span<T>& x) {
  return x.begin();
}

template <typename T>
auto end(Span<T>& x) {
  return x.end();
}

template <typename T>
auto end(const Span<T>& x) {
  return x.end();
}
} // namespace span
} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_SPAN_H
