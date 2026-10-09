#include "testDistrArray.h"

#include <molpro/linalg/array/DistrArrayGA.h>

using molpro::linalg::array::DistrArrayGA;

using ArrayTypes = ::testing::Types<DistrArrayGA>;
INSTANTIATE_TYPED_TEST_SUITE_P(GA, DistArrayBasicF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(GA, DistArrayBasicRMAF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(GA, DistrArrayRangeRMAF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(GA, DistrArrayRangeLinAlgF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(GA, DistrArrayRangeMinMaxF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(GA, TestDistrArray, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(GA, DistrArrayCollectiveLinAlgF, ArrayTypes);
