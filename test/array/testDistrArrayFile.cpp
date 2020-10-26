#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "parallel_util.h"

#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/Distribution.h>

using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::ScopeLock;
using molpro::linalg::test::mpi_comm;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

class DistrArrayFile_SetUp : public ::testing::Test {
  DistrArrayFile_SetUp() = default;
  void SetUp() override{};
  void TearDown() override{};
};

TEST(DistrArrayFile, constructor_default) {
  auto a = DistrArrayFile();
  ASSERT_TRUE(a.empty());
}

TEST(DistrArrayFile, constructor_dummy_with_filename) {
  auto a = DistrArrayFile(10, mpi_comm);
  EXPECT_EQ(a.size(), 10);
  ASSERT_FALSE(a.empty());
}

TEST(DistrArrayFile, constructor_fname_size) {
  {
    auto a = DistrArrayFile(100, mpi_comm);
    LockMPI3 lock{mpi_comm};
    {
      auto l = lock.scope();
      EXPECT_EQ(a.size(), 100);
      ASSERT_FALSE(a.empty());
    }
  }
}

TEST(DistrArrayFile, writeread) {
  ScopeLock l{mpi_comm};
  constexpr int n = 10;
  std::vector<double> v(n);
  std::iota(v.begin(), v.end(), 0.0);
  auto a = DistrArrayFile(10, mpi_comm);
  int mpi_size, mpi_rank;
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  MPI_Comm_size(mpi_comm, &mpi_size);
  auto dist = a.distribution();
  a.put(dist.range(mpi_rank).first, dist.range(mpi_rank).second, v.data());
  auto w{v};
  a.get(dist.range(mpi_rank).first, dist.range(mpi_rank).second, w.data());
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST(DistrArrayFile, copy) {
  ScopeLock l{mpi_comm};
  constexpr int n = 10;
  std::vector<double> v(n);
  std::iota(v.begin(), v.end(), 0.0);
  auto a = DistrArrayFile(10, mpi_comm);
  int mpi_size, mpi_rank;
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  MPI_Comm_size(mpi_comm, &mpi_size);
  auto dist = a.distribution();
  a.put(dist.range(mpi_rank).first, dist.range(mpi_rank).second, v.data());
  std::vector<double> w{v};
//  auto b = a; // TODO copy constructor not working yet
//  b.get(dist.range(mpi_rank).first, dist.range(mpi_rank).second, v.data());
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
//  b = a; // TODO operator=() not working yet
//  b.get(dist.range(mpi_rank).first, dist.range(mpi_rank).second, v.data());
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST(DistrArrayFile, move) {
  ScopeLock l{mpi_comm};
  constexpr int n = 10;
  std::vector<double> v(n);
  std::iota(v.begin(), v.end(), 0.0);
  auto a = DistrArrayFile(10, mpi_comm);
  int mpi_size, mpi_rank;
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  MPI_Comm_size(mpi_comm, &mpi_size);
  auto dist = a.distribution();
  a.put(dist.range(mpi_rank).first, dist.range(mpi_rank).second, v.data());
  DistrArrayFile b = std::move(a);
  auto w{v};
  b.get(dist.range(mpi_rank).first, dist.range(mpi_rank).second, w.data());
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST(DistrArrayFile, constructor_move) {
  auto&& a = DistrArrayFile{100, mpi_comm};
  {
    ScopeLock l{mpi_comm};
    DistrArrayFile b{std::move(a)};
    EXPECT_EQ(b.size(), 100);
    ASSERT_FALSE(b.empty());
  }
}

TEST(DistrArrayFile, DISABLED_compatible) {
  // TODO: CTEST randomly for parallel test
  auto a = DistrArrayFile{100, mpi_comm};
  auto b = DistrArrayFile{1000, mpi_comm};
  auto c = DistrArrayFile{};
  // auto d = DistrArrayFile{"test3.txt", 100, mpi_comm};
  ScopeLock l{mpi_comm};
  EXPECT_TRUE(a.compatible(a));
  EXPECT_TRUE(b.compatible(b));
  EXPECT_TRUE(c.compatible(c));
  // EXPECT_TRUE(a.compatible(d));
  EXPECT_FALSE(a.compatible(b));
  EXPECT_FALSE(a.compatible(c));
  EXPECT_FALSE(b.compatible(c));
  EXPECT_EQ(a.compatible(b), b.compatible(a));
  EXPECT_EQ(b.compatible(c), c.compatible(b));
  EXPECT_EQ(a.compatible(c), c.compatible(a));
  // EXPECT_EQ(a.compatible(d), d.compatible(a));
}
