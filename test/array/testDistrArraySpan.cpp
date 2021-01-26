#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "parallel_util.h"

#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/Distribution.h>

#ifdef LINEARALGEBRA_ARRAY_MPI3
#include <molpro/linalg/array/DistrArrayMPI3.h>
#endif

using molpro::linalg::array::Span;
using molpro::linalg::array::DistrArraySpan;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::ScopeLock;
using molpro::linalg::test::mpi_comm;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

class DistrArraySpan_Fixture : public ::testing::Test {
public:
  DistrArraySpan_Fixture() = default;
  void SetUp() override{
    a = DistrArraySpan(size, mpi_comm);
    auto dist = a.distribution();
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm, &mpi_size);
    left = dist.range(mpi_rank).first;
    right = dist.range(mpi_rank).second;
    for (int i = 0; i < mpi_size; i++) {
      displs.push_back(dist.range(i).first);
      chunks.push_back(dist.range(i).second - dist.range(i).first);
    }
  };
  void TearDown() override{};
  const size_t size = 12;
  int mpi_size, mpi_rank;
  int left, right;
  std::vector<int> chunks, displs;
  DistrArraySpan a;
};

TEST(DistrArraySpan, constructor_default) {
  auto a = DistrArraySpan();
  ASSERT_TRUE(a.empty());
}

TEST(DistrArraySpan, constructor_size) {
  {
    auto a = DistrArraySpan(100, mpi_comm);
    LockMPI3 lock{mpi_comm};
    {
      auto l = lock.scope();
      EXPECT_EQ(a.size(), 100);
      ASSERT_TRUE(a.empty()); // Different from DistrArrayMPI3 - OK?
    }
  }
}

TEST_F(DistrArraySpan_Fixture, constructor_copy) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  std::vector<double> w(size, 0);
  auto b = a;
  b.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArraySpan_Fixture, copy_assignment_empty) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  std::vector<double> w(size, 0);
  auto b = DistrArraySpan(size, mpi_comm);
  b = a;
  b.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArraySpan_Fixture, copy_assignment_filled) {
  std::vector<double> v(size), w(size);
  std::iota(v.begin(), v.end(), 0.5);
  std::iota(w.begin(), w.end(), 0.1);
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  auto b = DistrArraySpan(size, mpi_comm);
  b.allocate_buffer(Span<double>(&(*(w.begin() + left)), chunks[mpi_rank]));
  std::vector<double> x(size, 0);
  b = a;
  b.get(left, right, &(*(x.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, x.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(v, Pointwise(DoubleEq(), x));
}

TEST_F(DistrArraySpan_Fixture, span) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  std::vector<double> w(size, 0);
  v[0] += 0.5;
  a.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

#ifdef LINEARALGEBRA_ARRAY_MPI3
TEST(DistrArraySpan, constructor_copy_from_distr_array) {
  const double val = 0.5;
  auto a_mem = molpro::linalg::array::DistrArrayMPI3(100, mpi_comm);
  a_mem.allocate_buffer();
  a_mem.fill(val);
  auto a_disk = DistrArraySpan{a_mem};
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
#endif

TEST_F(DistrArraySpan_Fixture, constructor_move) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  DistrArraySpan b = std::move(a);
  std::vector<double> w(size, 0);
  b.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_EQ(b.size(), size);
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArraySpan_Fixture, assignment_move) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  auto b = DistrArraySpan();
  b = std::move(a);
  std::vector<double> w(size, 0);
  b.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_EQ(b.size(), size);
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST(DistrArraySpan, compatible) {
  auto a = DistrArraySpan{100, mpi_comm};
  auto b = DistrArraySpan{1000, mpi_comm};
  auto c = DistrArraySpan{};
  ScopeLock l{mpi_comm};
  EXPECT_TRUE(a.compatible(a));
  EXPECT_TRUE(b.compatible(b));
  EXPECT_TRUE(c.compatible(c));
  EXPECT_FALSE(a.compatible(b));
  EXPECT_FALSE(a.compatible(c));
  EXPECT_FALSE(b.compatible(c));
  EXPECT_EQ(a.compatible(b), b.compatible(a));
  EXPECT_EQ(b.compatible(c), c.compatible(b));
  EXPECT_EQ(a.compatible(c), c.compatible(a));
}

TEST_F(DistrArraySpan_Fixture, writeread) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  std::vector<double> w(size, 0);
  a.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArraySpan_Fixture, accumulate) {
  std::vector<double> v(size), w(size), x(size), y(size);
  std::iota(v.begin(), v.end(), 0.5);
  std::iota(w.begin(), w.end(), 0.5);
  std::transform(v.begin(), v.end(), w.begin(), x.begin(), [](auto& l, auto& r){return l + r;});
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  a.acc(left, right, &(*(w.begin() + left)));
  a.get(left, right, &(*(y.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(y, Pointwise(DoubleEq(), x));
}

TEST_F(DistrArraySpan_Fixture, gather) {
  std::vector<double> v(size), w(size);
  int n = -2;
  std::generate(v.begin(), v.end(), [&n]{ return n+=2; });
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  std::vector<DistrArraySpan::index_type> x(size/mpi_size);
  std::iota(x.begin(), x.end(), left);
  auto tmp = a.gather(x);
  for (int i = 0; i < x.size(); i++) {
    w[left + i] = tmp[i];
  }
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(w, Pointwise(DoubleEq(), v));
}

TEST_F(DistrArraySpan_Fixture, scatter) {
  std::vector<double> v(size, 0), w(size), y(size);
  int n = -2;
  std::generate(w.begin(), w.end(), [&n]{ return n+=2; });
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  std::vector<DistrArraySpan::index_type> x(size/mpi_size);
  std::iota(x.begin(), x.end(), left);
  std::vector<double> tmp(size/mpi_size);
  for (int i = 0; i < x.size(); i++) {
    tmp[i] = w[i + left];
  }
  a.scatter(x, tmp);
  a.get(left, right, &(*(y.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(y, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArraySpan_Fixture, scatter_acc) {
  std::vector<double> v(size), w(size), y(size);
  std::iota(v.begin(), v.end(), 0);
  int n = -2;
  std::generate(w.begin(), w.end(), [&n]{ return n+=2; });
  a.allocate_buffer(Span<double>(&(*(v.begin() + left)), chunks[mpi_rank]));
  std::vector<DistrArraySpan::index_type> x(size/mpi_size);
  std::iota(x.begin(), x.end(), left);
  std::vector<double> tmp(size/mpi_size);
  for (int i = 0; i < x.size(); i++) {
    tmp[i] = v[i + left];
  }
  a.scatter_acc(x, tmp);
  a.get(left, right, &(*(y.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(y, Pointwise(DoubleEq(), w));
}
