#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "parallel_util.h"

#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/Distribution.h>

#ifdef LINEARALGEBRA_ARRAY_MPI3
#include <molpro/linalg/array/DistrArrayMPI3.h>
#endif

using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::ScopeLock;
using molpro::linalg::test::mpi_comm;

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
  const size_t size = 12;
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
  std::vector<double> w(size, 0);
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
  for (int i = 0; i < x.size(); i++) {
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
  for (int i = 0; i < x.size(); i++) {
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
  for (int i = 0; i < x.size(); i++) {
    tmp[i] = v[i + left];
  }
  a.scatter_acc(x, tmp);
  a.get(left, right, &(*(y.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(y, Pointwise(DoubleEq(), w));
}
