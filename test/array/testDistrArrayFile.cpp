#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "parallel_util.h"

#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/array/default_handler.h>

#ifdef LINEARALGEBRA_ARRAY_MPI3
#include <molpro/linalg/array/DistrArrayMPI3.h>
using molpro::linalg::array::DistrArrayMPI3;
#endif

using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::DistrArraySpan;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::ScopeLock;
using molpro::linalg::test::mpi_comm;
using molpro::linalg::array::default_handler;
using molpro::linalg::array::ArrayHandlerDDiskDistr;
using molpro::linalg::array::Span;
using molpro::linalg::array::util::make_distribution_spread_remainder;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

class DistrArrayFile_Fixture : public ::testing::Test {
public:
  DistrArrayFile_Fixture() : a(DistrArrayFile(size, mpi_comm)) {}
  void SetUp() override {
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
  const size_t size = 1200;
  int mpi_size, mpi_rank;
  int left, right;
  std::vector<int> chunks, displs;
  DistrArrayFile a;
};

TEST(DistrArrayFile, constructor_size) {
  {
    auto a = DistrArrayFile(100, mpi_comm);
    LockMPI3 lock{mpi_comm};
    {
      auto l = lock.scope();
      EXPECT_EQ(a.size(), 100);
    }
  }
}

TEST_F(DistrArrayFile_Fixture, constructor_copy) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.put(left, right, &(*(v.cbegin() + left)));
  std::vector<double> w(size, 0), x(size, 0);
  auto b = a;
  b.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
  //  b = a; // TODO operator=() not working yet
  //  b.get(dist.range(mpi_rank).first, dist.range(mpi_rank).second, v.data());
  //  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

#ifdef LINEARALGEBRA_ARRAY_MPI3
TEST(DistrArrayFile, constructor_copy_from_distr_array) {
  const double val = 0.5;
  auto a_mem = molpro::linalg::array::DistrArrayMPI3(100, mpi_comm);
  a_mem.fill(val);
  auto a_disk = DistrArrayFile{a_mem};
  LockMPI3 lock{mpi_comm};
  auto vec = a_disk.vec();
  EXPECT_THAT(vec, Each(DoubleEq(val)));
  {
    auto l = lock.scope();
    EXPECT_EQ(a_disk.communicator(), a_mem.communicator());
    EXPECT_EQ(a_disk.size(), a_mem.size());
    EXPECT_TRUE(a_disk.distribution().compatible(a_mem.distribution()));
  }
}
#endif

TEST_F(DistrArrayFile_Fixture, handler_copies) {
  auto handler = ArrayHandlerDDiskDistr<DistrArrayFile, DistrArraySpan>{};
  size_t n = 3, dim = size;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim));
  std::vector<DistrArraySpan> mem_vecs;
  std::vector<DistrArrayFile> file_vecs;
  mem_vecs.reserve(n);
  file_vecs.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    mem_vecs.emplace_back(dim, Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength), comm_global());
    file_vecs.emplace_back(handler.copy(mem_vecs.back()));
  }
  for (size_t i = 0; i < n; i++) {
    file_vecs[i].get(left, right, &(*(vy[i].begin() + left)));
    //vy[i] = file_vecs[i].vec();
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, vy[i].data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
    EXPECT_THAT(vx[i], Pointwise(DoubleEq(), vy[i]));
  }
}

TEST_F(DistrArrayFile_Fixture, constructor_move) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.put(left, right, &(*(v.cbegin() + left)));
  DistrArrayFile b = std::move(a);
  std::vector<double> w(size, 0);
  b.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_EQ(b.size(), size);
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArrayFile_Fixture, assignment_move) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.put(left, right, &(*(v.cbegin() + left)));
  auto b = DistrArrayFile(1);
  b = std::move(a);
  std::vector<double> w(size, 0);
  b.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_EQ(b.size(), size);
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST(DistrArrayFile, compatible) {
  auto a = DistrArrayFile{100, mpi_comm};
  auto b = DistrArrayFile{1000, mpi_comm};
  ScopeLock l{mpi_comm};
  EXPECT_TRUE(a.compatible(a));
  EXPECT_TRUE(b.compatible(b));
  EXPECT_FALSE(a.compatible(b));
  EXPECT_EQ(a.compatible(b), b.compatible(a));
}

TEST_F(DistrArrayFile_Fixture, writeread) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.put(left, right, &(*(v.cbegin() + left)));
  std::vector<double> w(size, 0);
  a.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArrayFile_Fixture, accumulate) {
  std::vector<double> v(size), w(size), x(size), y(size);
  std::iota(v.begin(), v.end(), 0.5);
  std::iota(w.begin(), w.end(), 0.5);
  std::transform(v.begin(), v.end(), w.begin(), x.begin(), [](auto& l, auto& r) { return l + r; });
  a.put(left, right, &(*(v.cbegin() + left)));
  a.acc(left, right, &(*(w.begin() + left)));
  a.get(left, right, &(*(y.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(y, Pointwise(DoubleEq(), x));
}

TEST_F(DistrArrayFile_Fixture, gather) {
  std::vector<double> v(size), w(size);
  int n = -2;
  std::generate(v.begin(), v.end(), [&n] { return n += 2; });
  a.put(left, right, &(*(v.cbegin() + left)));
  std::vector<DistrArrayFile::index_type> x(size / mpi_size);
  std::iota(x.begin(), x.end(), left);
  auto tmp = a.gather(x);
  for (size_t i = 0; i < x.size(); i++) {
    w[left + i] = tmp[i];
  }
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(w, Pointwise(DoubleEq(), v));
}

TEST_F(DistrArrayFile_Fixture, scatter) {
  std::vector<double> v(size, 0), w(size), y(size);
  int n = -2;
  std::generate(w.begin(), w.end(), [&n] { return n += 2; });
  a.put(left, right, &(*(v.cbegin() + left)));
  std::vector<DistrArrayFile::index_type> x(size / mpi_size);
  std::iota(x.begin(), x.end(), left);
  std::vector<double> tmp(size / mpi_size);
  for (size_t i = 0; i < x.size(); i++) {
    tmp[i] = w[i + left];
  }
  a.scatter(x, tmp);
  a.get(left, right, &(*(y.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(y, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArrayFile_Fixture, scatter_acc) {
  std::vector<double> v(size), w(size), y(size);
  std::iota(v.begin(), v.end(), 0);
  int n = -2;
  std::generate(w.begin(), w.end(), [&n] { return n += 2; });
  a.put(left, right, &(*(v.cbegin() + left)));
  std::vector<DistrArrayFile::index_type> x(size / mpi_size);
  std::iota(x.begin(), x.end(), left);
  std::vector<double> tmp(size / mpi_size);
  for (size_t i = 0; i < x.size(); i++) {
    tmp[i] = v[i + left];
  }
  a.scatter_acc(x, tmp);
  a.get(left, right, &(*(y.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(y, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArrayFile_Fixture, dot_DistrArray) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0);
  a.put(left, right, &(*(v.cbegin() + left)));
  const DistrArraySpan s(size,Span<double>(&(*(v.begin() + left)),right-left));
  auto ss = s.dot(s);
  auto as = a.dot(s);
//  auto aa = a.dot(a);
  auto sa = s.dot(a);
  ScopeLock l{mpi_comm};
  EXPECT_NEAR(ss,size*(size-1)*(2*size-1)/6,1e-13);
  EXPECT_NEAR(as,size*(size-1)*(2*size-1)/6,1e-13);
//  EXPECT_NEAR(aa,size*(size-1)*(2*size-1)/6,1e-13);
  EXPECT_NEAR(sa,size*(size-1)*(2*size-1)/6,1e-13);
}
