#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "parallel_util.h"
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/DistrFlags.h>

using molpro::linalg::array::util::DistrFlags;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::test::mpi_comm;

namespace {
struct ScopeLock {
  ScopeLock() : lock{mpi_comm}, l{lock.scope()} {}

  LockMPI3 lock;
  decltype(std::declval<LockMPI3>().scope()) l;
};

int mpi_size() {
  int size;
  MPI_Comm_size(mpi_comm, &size);
  return size;
}
int mpi_rank() {
  int rank;
  MPI_Comm_rank(mpi_comm, &rank);
  return rank;
}
} // namespace

TEST(DistrFlags, constructor_comm) {
  DistrFlags df{mpi_comm};
  ScopeLock l{};
  ASSERT_FALSE(df.empty());
}

TEST(DistrFlags, constructor_comm_value) {
  DistrFlags df{mpi_comm, 0};
  ScopeLock l{};
  ASSERT_FALSE(df.empty());
}

TEST(DistrFlags, constructor_default) {
  DistrFlags df{};
  ScopeLock l{};
  ASSERT_TRUE(df.empty());
}
TEST(DistrFlags, constructor_copy) { DistrFlags df{}; }
TEST(DistrFlags, constructor_move) { DistrFlags df{}; }
TEST(DistrFlags, copy_assignment) { DistrFlags df{}; }
TEST(DistrFlags, move_assignment) { DistrFlags df{}; }

TEST(DistrFlags, access) {
  DistrFlags df{mpi_comm};
  auto p = df.access();
}
TEST(DistrFlags, access_rank) {
  auto rank = mpi_rank();
  auto size = mpi_size();
  DistrFlags df{mpi_comm};
  auto r = rank != size - 1 ? rank + 1 : 0;
//  ScopeLock l{};
  auto p = df.access(r);
}

// struct DistrFlagsF : DistrFlags, ::testing::Test {
//  DistrFlagsF() : DistrFlags(mpi_comm, value) {}
//
//  static constexpr int value = 1;
//};
