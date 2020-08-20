#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DEFAULT_HANDLER_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DEFAULT_HANDLER_H
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/ArrayHandlerDistr.h>
#include <molpro/linalg/array/ArrayHandlerDistrSparse.h>
#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/ArrayHandlerIterableSparse.h>
#include <molpro/linalg/array/ArrayHandlerSparse.h>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

template <typename T, typename = void>
struct is_iterable : std::false_type {};

template <typename T>
struct is_iterable<T, void_t<decltype(std::begin(std::declval<T>())), decltype(std::end(std::declval<T>()))>>
    : std::true_type {};

template <typename T, typename = void>
struct has_distr_tag : std::false_type {};

template <typename T>
struct has_distr_tag<T, void_t<typename T::distributed_array>> : std::true_type {};

} // namespace util

template <typename S, typename T, typename = std::enable_if_t<array::util::is_iterable<S>{}>,
          typename = std::enable_if_t<!array::util::has_mapped_type<S>{}>,
          typename = std::enable_if_t<array::util::is_iterable<T>{}>,
          typename = std::enable_if_t<!array::util::has_mapped_type<T>{}>>
auto create_default_handler() {
  return std::make_shared<array::ArrayHandlerIterable<S, T>>();
}

template <typename S, typename T, typename = std::enable_if_t<array::util::is_iterable<S>{}>,
          typename = std::enable_if_t<!array::util::has_mapped_type<S>{}>,
          typename = std::enable_if_t<array::util::has_mapped_type<T>{}>>
auto create_default_handler() {
  return std::make_shared<array::ArrayHandlerIterableSparse<S, T>>();
}

template <typename S, typename T, typename = std::enable_if_t<array::util::has_mapped_type<S>{}>,
          typename = std::enable_if_t<array::util::has_mapped_type<T>{}>>
auto create_default_handler() {
  return std::make_shared<array::ArrayHandlerSparse<S, T>>();
}

template <typename S, typename T, std::enable_if_t<array::util::has_distr_tag<S>{}, int> = 0,
          std::enable_if_t<array::util::has_distr_tag<T>{}, int> = 0>
auto create_default_handler() {
  return std::make_shared<array::ArrayHandlerDistr<S, T>>();
}

template <typename S, typename T, std::enable_if_t<array::util::has_distr_tag<S>{}, int> = 0,
          typename = std::enable_if_t<array::util::has_mapped_type<T>{}>>
auto create_default_handler() {
  return std::make_shared<array::ArrayHandlerDistrSparse<S, T>>();
}

} // namespace array
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DEFAULT_HANDLER_H
