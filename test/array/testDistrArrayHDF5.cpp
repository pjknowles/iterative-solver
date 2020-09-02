#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>

#include "data_util.h"
#include "parallel_util.h"

#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/PHDF5Handle.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/Distribution.h>

#ifdef LINEARALGEBRA_ARRAY_MPI3
#include <molpro/linalg/array/DistrArrayMPI3.h>
#endif

using molpro::linalg::array::DistrArrayHDF5;
using molpro::linalg::array::util::file_exists;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::PHDF5Handle;
using molpro::linalg::array::util::ScopeLock;
using molpro::linalg::test::mpi_comm;
using molpro::linalg::test::test_file_hdf5_n1;
using molpro::linalg::test::test_file_hdf5_n2;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

class DistrArrayHDF5_SetUp : public ::testing::Test {
public:
  DistrArrayHDF5_SetUp() = default;
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

TEST(DistrArrayHDF5, constructor_default) {
  auto a = DistrArrayHDF5{};
  ScopeLock l{mpi_comm};
  ASSERT_TRUE(a.empty());
}

TEST_F(DistrArrayHDF5_SetUp, constructor_dummy_with_fhandle) {
  auto a = DistrArrayHDF5{fhandle_n1};
  EXPECT_EQ(a.communicator(), fhandle_n1->communicator());
  EXPECT_EQ(a.size(), 0);
  EXPECT_TRUE(a.empty());
  EXPECT_FALSE(a.dataset_is_open());
}

TEST_F(DistrArrayHDF5_SetUp, constructor_fhandle_size) {
  ASSERT_TRUE(fhandle_n1);
  auto a = DistrArrayHDF5{fhandle_n1, size};
  LockMPI3 lock{mpi_comm};
  {
    auto l = lock.scope();
    EXPECT_EQ(a.communicator(), fhandle_n1->communicator());
    EXPECT_EQ(a.size(), size);
    EXPECT_EQ(a.file_handle(), fhandle_n1);
    EXPECT_TRUE(a.empty());
    EXPECT_FALSE(a.dataset_is_open());
  }
  a.close_access();
  {
    auto l = lock.scope();
    EXPECT_FALSE(a.dataset_is_open()) << "closing a closed file";
  }
  a.open_access();
  {
    auto l = lock.scope();
    EXPECT_TRUE(a.dataset_is_open());
  }
  a.close_access();
  {
    auto l = lock.scope();
    EXPECT_FALSE(a.dataset_is_open());
  }
}

TEST_F(DistrArrayHDF5_SetUp, constructor_move) {
  auto&& a = DistrArrayHDF5{fhandle_n1, size};
  DistrArrayHDF5 b{std::move(a)};
  ScopeLock l{mpi_comm};
  EXPECT_EQ(b.file_handle(), fhandle_n1);
  EXPECT_EQ(b.communicator(), fhandle_n1->communicator());
  EXPECT_EQ(b.size(), size);
  EXPECT_TRUE(b.empty());
}

TEST_F(DistrArrayHDF5_SetUp, compatible) {
  auto a = DistrArrayHDF5{fhandle_n1, size};
  auto b = DistrArrayHDF5{fhandle_n1};
  auto c = DistrArrayHDF5{};
  ScopeLock l{mpi_comm};
  EXPECT_TRUE(a.compatible(a));
  EXPECT_TRUE(b.compatible(b));
  EXPECT_TRUE(c.compatible(c));
  EXPECT_FALSE(a.compatible(b));
  EXPECT_FALSE(a.compatible(c));
  EXPECT_FALSE(a.compatible(c));
  EXPECT_FALSE(b.compatible(c));
  EXPECT_EQ(a.compatible(b), b.compatible(a));
  EXPECT_EQ(a.compatible(c), c.compatible(a));
}

#ifdef LINEARALGEBRA_ARRAY_MPI3
TEST_F(DistrArrayHDF5_SetUp, constructor_copy_from_distr_array) {
  const double val = 0.5;
  auto a_mem = molpro::linalg::array::DistrArrayMPI3(size, mpi_comm);
  a_mem.allocate_buffer();
  a_mem.fill(val);
  auto a_disk = DistrArrayHDF5{a_mem, fhandle_n1};
  LockMPI3 lock{mpi_comm};
  {
    auto l = lock.scope();
    EXPECT_EQ(a_disk.file_handle(), fhandle_n1);
    EXPECT_EQ(a_disk.communicator(), a_mem.communicator());
    EXPECT_EQ(a_disk.size(), a_mem.size());
    EXPECT_FALSE(a_disk.empty());
    EXPECT_TRUE(a_disk.distribution().compatible(a_mem.distribution()));
  }
  a_disk.open_access();
  {
    auto l = lock.scope();
    auto vec = a_disk.vec();
    EXPECT_THAT(vec, Each(DoubleEq(val)));
  }
}
#endif

TEST_F(DistrArrayHDF5_SetUp, CreateTempCopy) {
  LockMPI3 lock{mpi_comm};
  auto a = DistrArrayHDF5(fhandle_n1, size);
  std::string fname;
  {
    auto b = DistrArrayHDF5::CreateTempCopy(a);
    fname = b.file_handle()->file_name();
    b.open_access();
    auto l = lock.scope();
    ASSERT_TRUE(file_exists(fname));
  }
  auto l = lock.scope();
  ASSERT_FALSE(file_exists(fname));
}

struct DistrArrayHDF5_Fixture : DistrArrayHDF5_SetUp {
  void SetUp() override {
    DistrArrayHDF5_SetUp::SetUp();
    a = DistrArrayHDF5(fhandle_n1, size);
    a.open_access();
  }
  void TearDown() override {
    a.close_access();
    DistrArrayHDF5_SetUp::TearDown();
  }

  DistrArrayHDF5 a;
};

TEST_F(DistrArrayHDF5_Fixture, put_get) {
  auto ref_vec = std::vector<double>(size, 0.3);
  a.put(0, size, &ref_vec[0]);
  auto vec = a.get(0, size);
  ScopeLock l{mpi_comm};
  ASSERT_EQ(vec.size(), size);
  EXPECT_THAT(vec, Pointwise(DoubleEq(), ref_vec));
}

TEST_F(DistrArrayHDF5_Fixture, allocate_buffer_flush) {
  const double val = 11.;
  const double zero = 0.;
  a.fill(zero);
  a.allocate_buffer();
  LockMPI3 lock{mpi_comm};
  {
    auto l = lock.scope();
    ASSERT_FALSE(a.empty());
    auto buffer = a.local_buffer();
    a.fill(val);
    EXPECT_THAT(*buffer, Each(DoubleEq(val)));
  }
  auto vec = a.vec();
  {
    auto l = lock.scope();
    EXPECT_THAT(vec, Each(DoubleEq(zero)));
  }
  a.sync();
  a.flush();
  a.sync();
  a.free_buffer();
  vec = a.vec();
  {
    auto l = lock.scope();
    EXPECT_THAT(vec, Each(DoubleEq(val)));
  }
}
