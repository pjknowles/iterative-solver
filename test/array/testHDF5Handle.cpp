#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/HDF5Handle.h>

using molpro::linalg::array::util::file_exists;
using molpro::linalg::array::util::hdf5_open_file;

TEST(HDF5util, file_exists) { EXPECT_TRUE(file_exists(ARRAY_DATA)); }

TEST(HDF5util, hdf5_open_file_read_only) {
  auto id = hdf5_open_file(std::string(ARRAY_DATA) + "/single_dataset.hdf5");
  ASSERT_GE(id, 0) << "value <0 means file failed to open";
  H5Fclose(id);
}
