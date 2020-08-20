#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TYPE_TRAITS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TYPE_TRAITS_H
#include <type_traits>

namespace molpro {
namespace linalg {
namespace array {
namespace util {
//! Utility for metaprogramming that maps any types to void
template <class... Ts>
using void_t = void;

//! checks that type name A::mapped_type exists
template <class A, class = void_t<>>
struct has_mapped_type : std::false_type {};

template <class A>
struct has_mapped_type<A, void_t<typename A::mapped_type>> : std::true_type {};

//! Stores A::mapped_type or A::value_type as member type value, with former taking priority if both exist.
template <class A, bool = has_mapped_type<A>()>
struct mapped_or_value_type {
  using value = typename A::value_type;
};

template <class A>
struct mapped_or_value_type<A, true> {
  using value = typename A::mapped_type;
};

//! Checks that class T can be iterated with std::begin and std::end
template <class T, typename = void>
struct is_iterable : std::false_type {};

template <class T>
struct is_iterable<T, void_t<decltype(std::begin(std::declval<T>())), decltype(std::end(std::declval<T>()))>>
    : std::true_type {};

//! Checks if class T has a tag marking it as a distributed array
template <class T, typename = void>
struct has_distr_tag : std::false_type {};

template <class T>
struct has_distr_tag<T, void_t<typename T::distributed_array>> : std::true_type {};

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TYPE_TRAITS_H
