#include "testDistrArray.h"

#include <molpro/linalg/array/DistrArrayGA.h>

using molpro::linalg::array::DistrArrayGA;

using ArrayTypes = ::testing::Types<DistrArrayGA>;
INSTANTIATE_TYPED_TEST_SUITE_P(GA, TestDistrArray, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(GA, DistArrayInitializationF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(GA, DistrArrayRangeF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(GA, DistrArrayCollectiveOpF, ArrayTypes);
