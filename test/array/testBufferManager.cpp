#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "parallel_util.h"

#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/default_handler.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/BufferManager.h>
#include <molpro/linalg/array/util/Distribution.h>

#ifdef LINEARALGEBRA_ARRAY_HDF5
#include <molpro/linalg/array/DistrArrayHDF5.h>
using molpro::linalg::array::DistrArrayHDF5;
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
#include <molpro/linalg/array/DistrArrayGA.h>
using molpro::linalg::array::DistrArrayGA;
#endif
#ifdef LINEARALGEBRA_ARRAY_MPI3
#include <molpro/linalg/array/DistrArrayMPI3.h>
using molpro::linalg::array::DistrArrayMPI3;
#endif

using molpro::linalg::array::ArrayHandlerDDiskDistr;
using molpro::linalg::array::default_handler;
using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::DistrArraySpan;
using molpro::linalg::array::Span;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::make_distribution_spread_remainder;
using molpro::linalg::array::util::ScopeLock;
using molpro::linalg::test::mpi_comm;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;
using ::testing::Types;

template <class T>
class TestBufferManager : public ::testing::Test {
public:
  TestBufferManager() : f1(T(size, mpi_comm)), f2(T(size, mpi_comm)) {}
  void SetUp() override {
    auto dist = f1.distribution();
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm, &mpi_size);
    local_lo = dist.range(mpi_rank).first;
    local_hi = dist.range(mpi_rank).second;
    local_size = local_hi - local_lo;
    f_cvecref.emplace_back(std::cref(f1));
    f_cvecref.emplace_back(std::cref(f2));
    f_vector.push_back(f1);
    f_vector.push_back(f2);
    values1.resize(this->local_size);
    std::iota(values1.begin(), values1.end(), double(this->local_lo));
    values2.resize(this->local_size);
    std::iota(values2.begin(), values2.end(), double(this->local_lo + 100000));
  };
  void TearDown() override{};
  const size_t size = 1200;
  int mpi_size, mpi_rank;
  size_t local_lo, local_hi, local_size;
  T f1, f2;
  std::vector<std::reference_wrapper<const T>> f_cvecref;
  std::vector<T> f_vector;
  std::vector<double> values1;
  std::vector<double> values2;
};

TYPED_TEST_SUITE_P(TestBufferManager);

TYPED_TEST_P(TestBufferManager, Constructors) {
  molpro::linalg::array::util::BufferManager manager1(this->f_cvecref);
  molpro::linalg::array::util::BufferManager manager2(this->f_vector);
  molpro::linalg::array::util::BufferManager manager3(this->f1);
}

TYPED_TEST_P(TestBufferManager, TwoBuffers) {
  auto range = this->f1.distribution().range(molpro::mpi::rank_global());
  //  std::cout << "rank="<<molpro::mpi::rank_global()<<", range="<<this->local_lo<<":"<<this->local_hi<<std::endl;
  this->f1.put(this->local_lo, this->local_hi, this->values1.data());
  std::vector<double> received_values1(this->local_size);
  this->f1.get(this->local_lo, this->local_hi, received_values1.data());
  EXPECT_THAT(received_values1, Pointwise(DoubleEq(), this->values1));
  this->f2.put(this->local_lo, this->local_hi, this->values2.data());
  std::vector<double> received_values2(this->local_size);

  for (const auto& number_of_buffers : std::vector<size_t>{1, 2}) {
    for (const auto& buffer_size :
         std::vector<size_t>{1, this->local_size - 1, this->local_size, this->local_size + 1, this->local_size * 2}) {
      molpro::linalg::array::util::BufferManager manager(this->f_cvecref, buffer_size, number_of_buffers);
      std::fill(received_values1.begin(), received_values1.end(), -999);
      std::fill(received_values2.begin(), received_values2.end(), -999);
      for (const auto& buf : manager) {
        for (size_t i = 0; i < manager.buffer_size(); ++i)
          received_values1[manager.buffer_offset() + i] = buf[i];
        for (size_t i = 0; i < manager.buffer_size(); ++i)
          received_values2[manager.buffer_offset() + i] = buf[i + manager.buffer_stride()];
        EXPECT_LE(manager.buffer_size(), manager.buffer_stride())
            << "buffer_size " << manager.buffer_size() << ", buffer_stride " << manager.buffer_stride()
            << ", buffer_offset " << manager.buffer_offset() << std::endl;
        EXPECT_LE(manager.buffer_size() + manager.buffer_offset(), this->f1.size())
            << "buffer_size " << manager.buffer_size() << ", buffer_stride " << manager.buffer_stride()
            << ", buffer_offset " << manager.buffer_offset() << std::endl;
      }
      EXPECT_THAT(received_values1, Pointwise(DoubleEq(), this->values1));
      EXPECT_THAT(received_values2, Pointwise(DoubleEq(), this->values2));
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(TestBufferManager, TwoBuffers, Constructors);
using Implementations = ::testing::Types<
//#ifdef LINEARALGEBRA_ARRAY_HDF5 // TODO implement when DistrArrayHDF5 can take the constructor pattern
//    DistrArrayHDF5,
//#endif
#ifdef LINEARALGEBRA_ARRAY_GA
    DistrArrayGA,
#endif
#ifdef LINEARALGEBRA_ARRAY_MPI3
    DistrArrayMPI3,
#endif
    DistrArrayFile>;

INSTANTIATE_TYPED_TEST_SUITE_P(Implementations, TestBufferManager, Implementations);
