#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>

#include "testDistrArray.h"
#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/PHDF5Handle.h>
#include <molpro/linalg/array/util.h>

using molpro::linalg::array::DistrArrayHDF5;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

struct DistrArrayHDF5F : DistrArrayHDF5 {
  using DistrArrayHDF5::DistrArrayHDF5;
  DistrArrayHDF5F(size_t dim, MPI_Comm comm)
      : DistrArrayHDF5(std::make_shared<molpro::linalg::array::util::PHDF5Handle>("test_hdf5_array.h5", "/", comm),
                       dim) {
    DistrArrayHDF5::open_access();
  }

  ~DistrArrayHDF5F() override {
    DistrArrayHDF5::erase();
    DistrArrayHDF5::close_access();
  }

  void allocate_buffer() override{};
};

// Implement tests
