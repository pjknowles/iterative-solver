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
  constexpr int n=10;
  std::vector<double> v(n);
  std::iota(v.begin(),v.end(),0.0);
  auto a = DistrArrayFile(10, mpi_comm);
  a.put(0,n,v.data());
  std::vector<double> w(n);
  a.get(0,n,w.data());
  EXPECT_THAT(v, Pointwise(::testing::Eq(), w));
}

TEST(DistrArrayFile, DISABLED_copy) {
  constexpr int n=10;
  std::vector<double> v(n);
  std::iota(v.begin(),v.end(),0.0);
  auto a = DistrArrayFile(10, mpi_comm);
//  auto b = a;
  a.put(0,n,v.data());
  std::vector<double> w(n);
//  b.get(0,n,w.data());
  EXPECT_THAT(v, Pointwise(::testing::Eq(), w));
}

TEST(DistrArrayFile, move) {
  constexpr int n=10;
  std::vector<double> v(n);
  std::iota(v.begin(),v.end(),0.0);
  auto a = DistrArrayFile(10, mpi_comm);
  a.put(0,n,v.data());
  DistrArrayFile b = std::move(a);
  std::vector<double> w(n);
  b.get(0,n,w.data());
  EXPECT_THAT(v, Pointwise(::testing::Eq(), w));
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
  auto a = DistrArrayFile{ 100, mpi_comm};
  auto b = DistrArrayFile{ 1000, mpi_comm};
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
