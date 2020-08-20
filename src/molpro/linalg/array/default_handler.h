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

template <class T, class S, util::ArrayFamily = util::array_family_v<T>, util::ArrayFamily = util::array_family_v<S>>
struct default_handler {};

template <class T, class S>
struct default_handler<T, S, util::ArrayFamily::Iterable, util::ArrayFamily::Iterable> {
  using value = ArrayHandlerIterable<T, S>;
};

template <class T, class S>
struct default_handler<T, S, util::ArrayFamily::Sparse, util::ArrayFamily::Sparse> {
  using value = ArrayHandlerSparse<T, S>;
};

template <class T, class S>
struct default_handler<T, S, util::ArrayFamily::Distributed, util::ArrayFamily::Distributed> {
  using value = ArrayHandlerDistr<T, S>;
};

template <class T, class S>
struct default_handler<T, S, util::ArrayFamily::Iterable, util::ArrayFamily::Sparse> {
  using value = ArrayHandlerIterableSparse<T, S>;
};

template <class T, class S>
struct default_handler<T, S, util::ArrayFamily::Distributed, util::ArrayFamily::Sparse> {
  using value = ArrayHandlerDistrSparse<T, S>;
};

template <class T, class S>
using default_handler_t = typename default_handler<T, S>::value;

template <class T, class S>
auto create_default_handler() {
  return std::make_shared<default_handler_t<T, S>>();
}

} // namespace array
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DEFAULT_HANDLER_H
