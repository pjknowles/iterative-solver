#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/reset_dspace.h>

using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::itsolv::subspace::xspace::Dimensions;
using ::testing::DoubleEq;
using ::testing::Pointwise;
