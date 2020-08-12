#ifndef LINEARALGEBRA_TEST_ARRAY_TESTDISTRARRAYDISK_H
#define LINEARALGEBRA_TEST_ARRAY_TESTDISTRARRAYDISK_H
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>

#include "parallel_util.h"
#include <molpro/linalg/array/DistrArrayDisk.h>
#include <molpro/linalg/array/util.h>

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

#endif // LINEARALGEBRA_TEST_ARRAY_TESTDISTRARRAYDISK_H
