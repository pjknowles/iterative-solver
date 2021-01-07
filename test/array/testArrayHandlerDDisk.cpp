#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/array/PHDF5Handle.h>
#include <molpro/linalg/array/default_handler.h>
#include <molpro/linalg/array/util.h>

#include "data_util.h"
#include "parallel_util.h"

using molpro::linalg::array::ArrayHandlerDDisk;
using molpro::linalg::array::ArrayHandlerDDiskDistr;
using molpro::linalg::array::ArrayHandlerDistrDDisk;
using molpro::linalg::array::default_handler;
using molpro::linalg::array::DistrArrayHDF5;
using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::DistrArrayMPI3;
using molpro::linalg::array::util::file_exists;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::PHDF5Handle;
using molpro::linalg::test::mpi_comm;
using molpro::linalg::test::test_file_hdf5_n1;
using molpro::linalg::test::test_file_hdf5_n2;

using ::testing::ContainerEq;
using ::testing::Each;
using ::testing::DoubleEq;

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

TEST(ArrayHandlerDDisk, default_handler) {
  auto h = default_handler<DistrArrayHDF5, DistrArrayHDF5>::value{};
  ASSERT_TRUE((std::is_same<decltype(h), ArrayHandlerDDisk<DistrArrayHDF5>>::value));
}

TEST(ArrayHandlerDDisk_file, default_handler) {
  auto h = default_handler<DistrArrayFile, DistrArrayFile>::value{};
  ASSERT_TRUE((std::is_same<decltype(h), ArrayHandlerDDisk<DistrArrayFile>>::value));
}

TEST(ArrayHandlerDDiskDistr, default_handler) {
  auto h = default_handler<DistrArrayHDF5, DistrArrayMPI3>::value{};
  ASSERT_TRUE((std::is_same<decltype(h), ArrayHandlerDDiskDistr<DistrArrayHDF5, DistrArrayMPI3>>::value));
}

TEST(ArrayHandlerDDiskDistr_file, default_handler) {
  auto h = default_handler<DistrArrayFile, DistrArrayMPI3>::value{};
  ASSERT_TRUE((std::is_same<decltype(h), ArrayHandlerDDiskDistr<DistrArrayFile, DistrArrayMPI3>>::value));
}

TEST(ArrayHandlerDistrDDisk, default_handler) {
  auto h = default_handler<DistrArrayMPI3, DistrArrayHDF5>::value{};
  ASSERT_TRUE((std::is_same<decltype(h), ArrayHandlerDistrDDisk<DistrArrayMPI3, DistrArrayHDF5>>::value));
}

TEST(ArrayHandlerDistrDDisk_file, default_handler) {
  auto h = default_handler<DistrArrayMPI3, DistrArrayFile>::value{};
  ASSERT_TRUE((std::is_same<decltype(h), ArrayHandlerDistrDDisk<DistrArrayMPI3, DistrArrayFile>>::value));
}

TEST(ArrayHandlerDDiskDistr_File, constructor_copy_from_distr_array) {
  const double val = 0.5;
  auto a_mem = DistrArrayMPI3(100, mpi_comm);
  a_mem.allocate_buffer();
  a_mem.fill(val);
  //auto a_disk = DistrArrayFile(100, mpi_comm);
  auto h = default_handler<DistrArrayFile, DistrArrayMPI3>::value{};
  auto a_disk = h.copy(a_mem);
  LockMPI3 lock{mpi_comm};
  auto vec = a_disk.vec();
  EXPECT_THAT(vec, Each(DoubleEq(val)));
  {
    auto l = lock.scope();
    EXPECT_EQ(a_disk.communicator(), a_mem.communicator());
    EXPECT_EQ(a_disk.size(), a_mem.size());
    EXPECT_FALSE(a_disk.empty());
    EXPECT_TRUE(a_disk.distribution().compatible(a_mem.distribution()));
  }
}
