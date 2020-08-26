#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data_util.h"
#include "file_util.h"
#include "parallel_util.h"

#include <molpro/linalg/array/PHDF5Handle.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/temp_phdf5_handle.h>

using molpro::linalg::array::util::file_exists;
using molpro::linalg::array::util::HDF5Handle;
using molpro::linalg::array::util::PHDF5Handle;
using molpro::linalg::array::util::temp_phdf5_handle;
using molpro::linalg::test::GarbageCollector;
using molpro::linalg::test::name_inner_group_dataset;

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

struct PHDF5HandleTestFile : public ::testing::Test {
  PHDF5HandleTestFile() : file_name{"test_new_file.hdf5"}, lock(mpi_comm) { remove_file(); }
  ~PHDF5HandleTestFile() { remove_file(); }
  void remove_file() {
    int rank;
    MPI_Comm_rank(mpi_comm, &rank);
    if (rank == 0)
      if (file_exists(file_name))
        std::remove(file_name.c_str());
  }
  const std::string file_name;
  LockMPI3 lock;
};

TEST_F(PHDF5HandleTestFile, erase_on_destroy) {
  {
    auto handle = PHDF5Handle(file_name, mpi_comm);
    handle.open_file(PHDF5Handle::Access::read_write);
    auto l = lock.scope();
    ASSERT_FALSE(handle.erase_file_on_destroy());
    ASSERT_TRUE(handle.set_erase_file_on_destroy(true));
    ASSERT_TRUE(file_exists(handle.file_name()));
  }
  auto l = lock.scope();
  ASSERT_FALSE(file_exists(file_name));
}

TEST(temp_phdf5_handle, lifetime) {
  LockMPI3 lock{mpi_comm};
  auto g = GarbageCollector{};
  {
    auto h1 = temp_phdf5_handle(".temp", mpi_comm);
    {
      auto l = lock.scope();
      ASSERT_FALSE(h1.file_name().empty());
      ASSERT_FALSE(file_exists(h1.file_name()));
      ASSERT_TRUE(h1.erase_file_on_destroy());
    }
    g.file_name = h1.file_name();
    h1.open_file(HDF5Handle::Access::read_write);
    {
      auto l = lock.scope();
      ASSERT_TRUE(file_exists(h1.file_name()));
    }
  }
  auto l = lock.scope();
  ASSERT_FALSE(file_exists(g.file_name));
}
