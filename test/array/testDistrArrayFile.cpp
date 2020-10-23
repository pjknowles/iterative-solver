#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>

#include "parallel_util.h"

#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/Distribution.h>

using molpro::linalg::test::mpi_comm;
using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::ScopeLock;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

class DistrArrayFile_SetUp : public ::testing::Test {
  DistrArrayFile_SetUp() = default;
  void SetUp() override {
  };
  void TearDown() override {};
};

TEST(DistrArrayFile, constructor_default) {
  auto a = DistrArrayFile();
  ASSERT_TRUE(a.empty());
}

TEST(DistrArrayFile, constructor_dummy_with_filename) {
  {
    auto a = DistrArrayFile("test.txt", mpi_comm);
    EXPECT_EQ(a.size(), 0);
    ASSERT_TRUE(a.empty());
  }
  ASSERT_FALSE(std::fstream{"test.txt"}.is_open());
  ASSERT_TRUE(std::fstream{"test.txt"}.fail());
}

TEST(DistrArrayFile, constructor_fname_size) {
  {
    auto a = DistrArrayFile("test.txt", 100, mpi_comm);
    LockMPI3 lock{mpi_comm};
    {
      auto l = lock.scope();
      EXPECT_EQ(a.size(), 100);
      ASSERT_TRUE(a.empty());
    }
    a.close_access();
    {
      auto l = lock.scope();
      ASSERT_FALSE(a.is_file_open());
    }
    a.open_access();
    {
      auto l = lock.scope();
      ASSERT_TRUE(a.is_file_open());
    }
    a.erase();
    {
      auto l = lock.scope();
      ASSERT_TRUE(a.is_file_open());
    }
  }
  ASSERT_FALSE(std::fstream{"test.txt"}.is_open());
  ASSERT_TRUE(std::fstream{"test.txt"}.fail());
}

TEST(DistrArrayFile, constructor_move) {
  auto&& a = DistrArrayFile{"test.txt", 100, mpi_comm};
  {
    ScopeLock l{mpi_comm};
    DistrArrayFile b{std::move(a)};
    EXPECT_EQ(b.size(), 100);
    ASSERT_TRUE(b.empty());
  }
  ASSERT_FALSE(std::fstream{"test.txt"}.is_open());
  ASSERT_TRUE(std::fstream{"test.txt"}.fail());
}

TEST(DistrArrayFile, DISABLED_compatible) {
  //TODO: CTEST randomly for parallel test
  auto a = DistrArrayFile{"test.txt", 100, mpi_comm};
  auto b = DistrArrayFile{"test2.txt", mpi_comm};
  auto c = DistrArrayFile{};
  //auto d = DistrArrayFile{"test3.txt", 100, mpi_comm};
  ScopeLock l{mpi_comm};
  EXPECT_TRUE(a.compatible(a));
  EXPECT_TRUE(b.compatible(b));
  EXPECT_TRUE(c.compatible(c));
  //EXPECT_TRUE(a.compatible(d));
  EXPECT_FALSE(a.compatible(b));
  EXPECT_FALSE(a.compatible(c));
  EXPECT_FALSE(b.compatible(c));
  EXPECT_EQ(a.compatible(b), b.compatible(a));
  EXPECT_EQ(b.compatible(c), c.compatible(b));
  EXPECT_EQ(a.compatible(c), c.compatible(a));
  //EXPECT_EQ(a.compatible(d), d.compatible(a));
}
