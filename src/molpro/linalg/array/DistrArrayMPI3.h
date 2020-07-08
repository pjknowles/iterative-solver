#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYMPI3_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYMPI3_H
#include "molpro/gci/array/DistrArray.h"

namespace molpro {
namespace gci {
namespace array {

/*!
 * @brief Implementation of distributed array using MPI3 RMA operations
 *
 * The array buffer is a window which is open for all processes on creation and
 * remains open until the destruction of the array.
 * @warning Care must be taken that overlapping put and get operations or linear algebra
 * do not cause undefined behaviour.
 */
class DistrArrayMPI3 : public DistrArray {};

} // namespace array
} // namespace gci
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYMPI3_H
