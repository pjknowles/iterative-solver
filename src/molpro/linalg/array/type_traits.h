#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TYPE_TRAITS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TYPE_TRAITS_H
#include <cmath>
#include <cstdlib>
#include <type_traits>

namespace molpro {
namespace linalg {
namespace array {
//! Utility for metaprogramming that maps any types to void
template <class... Ts>
using void_t = void;

//! checks that type name A::mapped_type exists
template <class A, class = void_t<>>
struct has_mapped_type : std::false_type {};

template <class A>
struct has_mapped_type<A, void_t<typename A::mapped_type>> : std::true_type {};

template <class A>
constexpr bool has_mapped_type_v = has_mapped_type<A>{};

template <class T>
constexpr bool is_sparse_v = has_mapped_type_v<T>;

//! Stores A::mapped_type or A::value_type as member type value, with former taking priority if both exist.
template <class A, bool = has_mapped_type_v<A>>
struct mapped_or_value_type {
  using type = typename A::value_type;
};

template <class A>
struct mapped_or_value_type<A, true> {
  using type = typename A::mapped_type;
};

template <class A>
using mapped_or_value_type_t = typename mapped_or_value_type<A>::type;

//! Checks that class T can be iterated with std::begin and std::end, and is not sparse
template <class T, typename = void>
struct is_iterable : std::false_type {};

template <class T>
struct is_iterable<T, void_t<decltype(std::begin(std::declval<T>())), decltype(std::end(std::declval<T>())),
                             std::enable_if_t<!is_sparse_v<T>>>> : std::true_type {};

template <class T>
constexpr bool is_iterable_v = is_iterable<T>{};

//! Checks if class T has a tag marking it as a distributed array
template <class T, typename = void>
struct is_distributed : std::false_type {};

template <class T>
struct is_distributed<T, void_t<typename T::distributed_array>> : std::true_type {};

template <class T>
constexpr bool is_distributed_v = is_distributed<T>::value;

//! A tag to distinguish different families for array types, e.g. std::vector<> is iterable, std::map is sparse etc.
enum class ArrayFamily { None, Iterable, Sparse, Distributed };

//! Deduces which family an array type belongs to
template <class T, bool = is_iterable_v<T>, bool = is_sparse_v<T>, bool = is_distributed_v<T>>
struct array_family {
  constexpr auto value() { return ArrayFamily::None; }
};

template <class T>
struct array_family<T, true, false, false> {
  constexpr auto value() { return ArrayFamily::Iterable; }
};

template <class T>
struct array_family<T, false, true, false> {
  constexpr auto value() { return ArrayFamily::Sparse; }
};

template <class T>
struct array_family<T, false, false, true> {
  constexpr auto value() { return ArrayFamily::Distributed; }
};

template <class T>
constexpr auto array_family_v = array_family<T>{}.value();

template <typename T>
constexpr auto check_abs() {
  using std::abs;
  return abs(T{});
}
} // namespace array
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TYPE_TRAITS_H