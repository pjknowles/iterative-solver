#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "parallel_util.h"
#include <molpro/linalg/array/util.h>

using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::test::mpi_comm;

TEST(LockMPI3, creation) { LockMPI3 lock{mpi_comm}; }

TEST(LockMPI3, lock_unlock) {
  LockMPI3 lock{mpi_comm};
  ASSERT_NO_FATAL_FAILURE(lock.lock());
  ASSERT_NO_FATAL_FAILURE(lock.unlock());
}

TEST(LockMPI3, scope) {
  LockMPI3 lock{mpi_comm};
  auto proxy = lock.scope();
}

TEST(LockMPI3, scope_and_delete) {
  auto lock = std::make_shared<LockMPI3>(mpi_comm);
  auto proxy = lock->scope();
  lock.reset();
}

// Simulate Fetch_and_op with put and get. If the lock mechanism does not work this is very likely to fail
// although it is not guranteed
TEST(LockMPI3, locking_mechanism_simulates_fetch_and_op) {
  auto comm = mpi_comm;
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  if (size == 1)
    return;
  LockMPI3 lock{comm};
  MPI_Win win = MPI_WIN_NULL;
  int *base = nullptr;
  MPI_Win_allocate(1*sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &base, &win);
  MPI_Win_fence(0, win);
  int zero = 0;
  MPI_Put(&zero, 1, MPI_INT, rank, 0, 1, MPI_INT, win);
  MPI_Win_fence(0, win);
  MPI_Barrier(comm);
  int v = -1;
  {
    auto proxy = lock.scope();
    MPI_Get(&v, 1, MPI_INT, 0, 0, 1, MPI_INT, win);
    ++v;
    MPI_Put(&v, 1, MPI_INT, 0, 0, 1, MPI_INT, win);
  }
  MPI_Win_fence(0, win);
  if (rank == 0)
    v = base[0];
  MPI_Bcast(&v, 1, MPI_INT, 0, comm);
  {
    auto proxy = lock.scope();
    ASSERT_EQ(v, size);
  }
  MPI_Win_free(&win);
}