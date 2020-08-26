#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "file_util.h"
#include <molpro/linalg/array/util/temp_file.h>
#include <molpro/linalg/array/util/temp_hdf5_handle.h>

using molpro::linalg::array::util::file_exists;
using molpro::linalg::array::util::HDF5Handle;
using molpro::linalg::array::util::temp_file_name;
using molpro::linalg::array::util::temp_hdf5_handle;
using molpro::linalg::test::GarbageCollector;

TEST(TempFile, temp_file_name) {
  const auto body = ".temp";
  const auto suffix = ".suffix";
  auto f1 = temp_file_name(body, suffix);
  auto f2 = temp_file_name(body, suffix);
  ASSERT_FALSE(f1.empty());
  ASSERT_FALSE(file_exists(f1));
  ASSERT_EQ(f1, f2);
  std::ofstream(f1.c_str()).close();
  auto g = GarbageCollector(f1);
  ASSERT_TRUE(file_exists(f1));
  auto f3 = temp_file_name(body, suffix);
  ASSERT_FALSE(f1.empty());
  ASSERT_NE(f1, f3);
}
