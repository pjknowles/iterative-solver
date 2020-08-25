#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/util/TempHandle.h>
#include <molpro/linalg/array/util/temp_file.h>
#include <molpro/linalg/array/util/temp_hdf5_handle.h>

using molpro::linalg::array::util::file_exists;
using molpro::linalg::array::util::temp_file_name;
using molpro::linalg::array::util::temp_hdf5_handle;

class GarbageCollector {
public:
  GarbageCollector(std::string fname) : file_name(std::move(fname)) { remove_test_files(); }
  ~GarbageCollector() { remove_test_files(); }
  void remove_test_files() {
    if (file_exists(file_name))
      std::remove(file_name.c_str());
  }

  std::string file_name;
};

TEST(TempFile, temp_file_name) {
  auto f1 = temp_file_name(".temp", "");
  ASSERT_FALSE(f1.empty());
  ASSERT_FALSE(file_exists(f1));
  auto f2 = temp_file_name(".temp", "");
  ASSERT_FALSE(f1.empty());
  ASSERT_FALSE(file_exists(f1));
  ASSERT_NE(f1, f2);
}

TEST(TempFile, temp_hdf5_handle) { auto h1 = temp_hdf5_handle(".temp"); }
