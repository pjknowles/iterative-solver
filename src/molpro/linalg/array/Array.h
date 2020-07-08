#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_ARRAY_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_ARRAY_H

#ifdef GCI_ARRAY_ArrayGA
#define GCI_ARRAY_ARRAY_TYPE ArrayGA
#include "molpro/gci/array/ArrayGA.h"
#elif GCI_ARRAY_ArrayMPI3
#define GCI_ARRAY_ARRAY_TYPE ArrayMPI3
#include "molpro/gci/array/ArrayMPI3.h"
#else
#define GCI_ARRAY_ARRAY_TYPE ArrayMPI3
#include "molpro/gci/array/ArrayMPI3.h"
#endif

namespace molpro {
namespace gci {
namespace array {
using Array = GCI_ARRAY_ARRAY_TYPE;
} // namespace array
} // namespace gci
} // namespace molpro
#undef GCI_ARRAY_ARRAY_TYPE
#endif // GCI_SRC_MOLPRO_GCI_ARRAY_ARRAY_H
