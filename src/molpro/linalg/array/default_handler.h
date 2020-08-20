#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DEFAULT_HANDLER_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DEFAULT_HANDLER_H
#include <molpro/linalg/array/ArrayHandlerDistr.h>
#include <molpro/linalg/array/ArrayHandlerDistrSparse.h>
#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/ArrayHandlerIterableSparse.h>
#include <molpro/linalg/array/ArrayHandlerSparse.h>
#include <molpro/linalg/array/util/type_traits.h>

namespace molpro {
namespace linalg {
namespace array {

template <typename S, typename T, typename = std::enable_if_t<util::is_iterable_v<S>>,
          typename = std::enable_if_t<!util::has_mapped_type_v<S>>, typename = std::enable_if_t<util::is_iterable_v<T>>,
          typename = std::enable_if_t<!util::has_mapped_type_v<T>>>
auto create_default_handler() {
  return std::make_shared<ArrayHandlerIterable<S, T>>();
}

template <typename S, typename T, typename = std::enable_if_t<util::is_iterable_v<S>>,
          typename = std::enable_if_t<!util::has_mapped_type_v<S>>,
          typename = std::enable_if_t<util::has_mapped_type_v<T>>>
auto create_default_handler() {
  return std::make_shared<ArrayHandlerIterableSparse<S, T>>();
}

template <typename S, typename T, typename = std::enable_if_t<util::has_mapped_type_v<S>>,
          typename = std::enable_if_t<util::has_mapped_type_v<T>>>
auto create_default_handler() {
  return std::make_shared<ArrayHandlerSparse<S, T>>();
}

template <typename S, typename T, std::enable_if_t<util::has_distr_tag_v<S>, int> = 0,
          std::enable_if_t<util::has_distr_tag_v<T>, int> = 0>
auto create_default_handler() {
  return std::make_shared<ArrayHandlerDistr<S, T>>();
}

template <typename S, typename T, std::enable_if_t<util::has_distr_tag_v<S>, int> = 0,
          typename = std::enable_if_t<util::has_mapped_type_v<T>>>
auto create_default_handler() {
  return std::make_shared<ArrayHandlerDistrSparse<S, T>>();
}

} // namespace array
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DEFAULT_HANDLER_H
