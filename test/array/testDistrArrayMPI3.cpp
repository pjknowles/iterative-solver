#include "testDistrArray.h"

#include <molpro/linalg/array/DistrArrayMPI3.h>

using molpro::linalg::array::DistrArrayMPI3;

using ArrayTypes = ::testing::Types<DistrArrayMPI3>;
INSTANTIATE_TYPED_TEST_SUITE_P(MPI3, TestDistrArray, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(MPI3, DistArrayInitializationF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(MPI3, DistrArrayRangeF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(MPI3, DistrArrayCollectiveOpF, ArrayTypes);