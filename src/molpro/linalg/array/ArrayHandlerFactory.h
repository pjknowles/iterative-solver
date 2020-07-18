#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERFACTORY_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERFACTORY_H
#include "molpro/linalg/array/ArrayHandler.h"
#include "molpro/linalg/array/ArrayHandlerIterable.h"
#include "molpro/linalg/array/DistrArray.h"
#include "molpro/linalg/array/DistrArrayHDF5.h"

namespace molpro {
namespace linalg {
namespace array {

/*!
 * @brief Factory for creating correct handler for any combination of LinearAlgebra arrays and iterable containers
 *
 *
 * I need to be able to generate correct ArrayHandler for any pair of arrays
 * We will have the following types of arrays:
 *   * Iterable arrays - any container that has begin/cbegin(), end()/cend(), push_back and copy constructor Array(const
 *                       Array&)
 *   * Disk arrays - array residing on disk in an HDF5 file
 *   * Distributed arrays in core - array distributed among processors in a PGAS model
 *   * Distributed arrays on disk - distributed array stored on disk
 *
 * ArrayHandler should be symmetric, order of arrays shouldn't matter.
 *
 * ArrayHandlers that we will provide:
 *   1. <Iterable, Iterable> - pair of iterable arrays, not necessarily contiguous
 *   2. <ContiguousIterable, Disk> - contiguous array that is iterable and has contiguous data buffer accessed by the
 * .data() method, paired with an array on disk (e.g. ArrayHDF5)
 *   3. <DistrArray, DistrArray> - pair of DistrArrays, implementation type does not matter
 *   4. <DistrArray, DistrArrayHDF5> - any DistrArray paired with an HDF5 version for dumping it to disk.
 *
 * Priority for ArrayHandler selection:
 *  - From most specialised to least specialised,
 *    4 > 3 > 2 > 1
 *    In a tree form,
 *         1 -> 2
 *           -> 3 -> 4
 *
 * Example usage
 * -------------
 * @code{.cpp}
 * auto array_handler = ArrayHandlerFactory<DistrArrayMPI3, DistrArrayHDF5>()::create();
 * @endcode
 *
 */
template <typename AL, typename AR = AL> struct ArrayHandlerFactory {
  //! Pathway when both arrays are derived from DistrArray
  template <typename T1 = AL, typename T2 = AR,
            typename = typename std::enable_if_t<std::is_base_of<DistrArray, T1>::value &&
                                                 std::is_base_of<DistrArray, T2>::value>>
  static auto create() {
    return createDistrArray<T1, T2>();
  }

  //! Pathway when neither arrays are derived from DistrArray
  template <
      typename T1 = AL, typename T2 = AR,
      typename std::enable_if_t<!std::is_base_of<DistrArray, T1>::value && !std::is_base_of<DistrArray, T2>::value,
                                nullptr_t> = nullptr>
  static auto create() {
    return createIterable<T1, T2>();
  }

  //! handler for a pair of Iterable arrays
  template <typename T1, typename T2> static auto createIterable() {
    return std::make_shared<ArrayHandlerIterable<T1, T2>>();
  }

  //! handler for pair of DistrArray
  template <typename T1, typename T2,
            typename = typename std::enable_if_t<!std::is_base_of<DistrArrayHDF5, T1>::value &&
                                                 !std::is_base_of<DistrArrayHDF5, AR>::value>>
  static auto createDistrArray() {
    return;
  }

  //! handler for <DistrArrayHDF5 , DistrArray> and  <DistrArray, DistrArrayHDF5>
  template <typename T1, typename T2,
            typename std::enable_if_t<
                std::is_base_of<DistrArrayHDF5, T1>::value || std::is_base_of<DistrArrayHDF5, AR>::value, nullptr_t> =
                nullptr>
  static auto createDistrArray() {
    return;
  }
};

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERFACTORY_H
