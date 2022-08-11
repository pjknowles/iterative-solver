#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "parallel_util.h"

#include <molpro/linalg/array/Resource.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/wrap.h>


using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

TEST(TestResource, construct) {
  molpro::linalg::array::Resource resource;
  EXPECT_EQ(resource.m_mpi_communicator , molpro::mpi::comm_global());
  EXPECT_EQ(resource.m_directory , ".");
  molpro::linalg::array::Resource resourc2{molpro::mpi::comm_self(),"/etc"};
  EXPECT_EQ(resourc2.m_mpi_communicator , molpro::mpi::comm_self());
  EXPECT_EQ(resourc2.m_directory , "/etc");
}