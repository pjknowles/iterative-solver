#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data_util.h"
#include "parallel_util.h"

#include <molpro/linalg/array/PHDF5Handle.h>
#include <molpro/linalg/array/util.h>

using molpro::linalg::array::util::HDF5Handle;
using molpro::linalg::array::util::PHDF5Handle;

using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::test::mpi_comm;

TEST(PHDF5Handle, constructor) {
  auto h = PHDF5Handle{mpi_comm};
  EXPECT_EQ(h.communicator(), mpi_comm);
}

TEST(PHDF5Handle, open_and_close_file) {
  LockMPI3 lock{mpi_comm};
  auto h = PHDF5Handle{name_inner_group_dataset, mpi_comm};
  auto fid = h.open_file(HDF5Handle::Access::read_only);
  {
    auto l = lock.scope();
    EXPECT_NE(fid, HDF5Handle::hid_default);
    EXPECT_TRUE(h.file_is_open());
  }
  h.close_file();
}

TEST(PHDF5Handle, open_and_close_group) {
  LockMPI3 lock{mpi_comm};
  auto group_name = "/group1";
  auto h = PHDF5Handle{name_inner_group_dataset, mpi_comm};
  auto gid = h.open_group(group_name);
  {
    auto l = lock.scope();
    EXPECT_EQ(gid, h.group_id());
    EXPECT_NE(h.group_id(), HDF5Handle::hid_default);
    EXPECT_NE(h.file_id(), HDF5Handle::hid_default);
    EXPECT_TRUE(h.file_is_open());
    EXPECT_TRUE(h.group_is_open());
  }
  h.close_group();
  {
    auto l = lock.scope();
    EXPECT_EQ(h.group_id(), HDF5Handle::hid_default);
    EXPECT_NE(h.file_id(), HDF5Handle::hid_default);
    EXPECT_TRUE(h.file_is_open());
    EXPECT_FALSE(h.group_is_open());
  }
  h.close_file();
  {
    auto l = lock.scope();
    EXPECT_EQ(h.group_id(), HDF5Handle::hid_default);
    EXPECT_EQ(h.file_id(), HDF5Handle::hid_default);
    EXPECT_FALSE(h.file_is_open());
    EXPECT_FALSE(h.group_is_open());
  }
}
