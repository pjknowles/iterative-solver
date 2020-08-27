#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/ArrayHandlerDDisk.h>
#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/PHDF5Handle.h>
#include <molpro/linalg/array/util.h>

#include "data_util.h"
#include "parallel_util.h"

using molpro::linalg::array::ArrayHandlerDDisk;
using molpro::linalg::array::DistrArrayHDF5;
using molpro::linalg::array::util::file_exists;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::PHDF5Handle;
using molpro::linalg::test::mpi_comm;
using molpro::linalg::test::test_file_hdf5_n1;
using molpro::linalg::test::test_file_hdf5_n2;

using ::testing::ContainerEq;
using ::testing::Each;

class ArrayHandlerDDiskF : public ::testing::Test {
public:
  ArrayHandlerDDiskF() = default;
  void SetUp() override {
    remove_test_files();
    fhandle_n1 = std::make_shared<PHDF5Handle>(test_file_hdf5_n1, "/", mpi_comm);
    fhandle_n2 = std::make_shared<PHDF5Handle>(test_file_hdf5_n2, "/", mpi_comm);
  }

  void TearDown() override { remove_test_files(); }
  static void remove_test_files() {
    int rank;
    MPI_Comm_rank(mpi_comm, &rank);
    if (rank == 0) {
      for (const auto& f : {test_file_hdf5_n1, test_file_hdf5_n2})
        if (file_exists(f))
          std::remove(f.c_str());
    }
    MPI_Barrier(mpi_comm);
  }

  std::shared_ptr<PHDF5Handle> fhandle_n1;
  std::shared_ptr<PHDF5Handle> fhandle_n2;
  const size_t size = 30;
};

TEST_F(ArrayHandlerDDiskF, constructor_default) {
  LockMPI3 lock{mpi_comm};
  auto a = ArrayHandlerDDisk<DistrArrayHDF5>();
  auto x = DistrArrayHDF5(fhandle_n1);
  auto y = a.copy(x);
  ASSERT_NE(y.file_handle(), x.file_handle());
  ASSERT_TRUE(y.file_handle()->erase_file_on_destroy());
}

TEST_F(ArrayHandlerDDiskF, constructor_erase_on_default) {
  LockMPI3 lock{mpi_comm};
  auto copy_func = [this](const auto& source) {
    fhandle_n2->set_erase_file_on_destroy(true);
    return DistrArrayHDF5(source, fhandle_n2);
  };
  auto a = ArrayHandlerDDisk<DistrArrayHDF5>(copy_func);
  auto x = DistrArrayHDF5(fhandle_n1);
  {
    auto y = a.copy(x);
    auto l = lock.scope();
    EXPECT_TRUE(fhandle_n2->erase_file_on_destroy());
    ASSERT_EQ(y.file_handle(), fhandle_n2);
  }
  {
    auto l = lock.scope();
    ASSERT_FALSE(file_exists(fhandle_n2->file_name()));
  }
}
