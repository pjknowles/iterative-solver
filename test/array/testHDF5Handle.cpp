#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/HDF5Handle.h>

using molpro::linalg::array::util::file_exists;
using molpro::linalg::array::util::hdf5_link_exists;
using molpro::linalg::array::util::hdf5_open_file;
using molpro::linalg::array::util::HDF5Handle;

namespace {
// File contents: "/dataset" where dataset is a float64 arange(0,30)
const std::string name_single_dataset{std::string(ARRAY_DATA) + "/single_dataset.hdf5"};
// File contents: "/group1/group2/dataset" where dataset is a float64 arange(0,30)
const std::string name_inner_group_dataset{std::string(ARRAY_DATA) + "/inner_group_dataset.hdf5"};
} // namespace

TEST(HDF5util, file_exists) { EXPECT_TRUE(file_exists(ARRAY_DATA)); }

TEST(HDF5util, hdf5_open_file_read_only) {
  auto id = hdf5_open_file(name_single_dataset);
  ASSERT_GE(id, 0) << "value <0 means file failed to open";
  H5Fclose(id);
}

class HDF5utilF : public ::testing::Test {
public:
  HDF5utilF() = default;
  void SetUp() override {
    id_single = hdf5_open_file(name_single_dataset);
    ASSERT_GE(id_single, 0) << "value <0 means file failed to open";
    id_nested = hdf5_open_file(name_inner_group_dataset);
    ASSERT_GE(id_nested, 0) << "value <0 means file failed to open";
  }
  void TearDown() override {
    H5Fclose(id_single);
    H5Fclose(id_nested);
  }

  hid_t id_single = -1;
  hid_t id_nested = -1;
};

TEST_F(HDF5utilF, hdf5_link_exists) {
  EXPECT_LE(hdf5_link_exists(id_single, ""), 0);
  EXPECT_GT(hdf5_link_exists(id_single, "/dataset"), 0);
  EXPECT_GT(hdf5_link_exists(id_single, " /dataset"), 0);
  EXPECT_GT(hdf5_link_exists(id_single, " /dataset/ "), 0);
  EXPECT_GT(hdf5_link_exists(id_nested, "/group1"), 0);
  EXPECT_GT(hdf5_link_exists(id_nested, "/group1/group2"), 0);
  EXPECT_GT(hdf5_link_exists(id_nested, "/group1/group2/dataset"), 0);
  EXPECT_LE(hdf5_link_exists(id_nested, ""), 0);
  EXPECT_LE(hdf5_link_exists(id_nested, "does_not_exist"), 0);
  EXPECT_LE(hdf5_link_exists(id_nested, "/does_not_exist"), 0);
  EXPECT_LE(hdf5_link_exists(id_nested, "/group1/does_not_exist"), 0);
}

TEST(HDF5HandleDefault, default_constructor) {
  auto h = HDF5Handle();
  ASSERT_TRUE(h.empty());
  EXPECT_FALSE(h.file_owner());
  EXPECT_FALSE(h.group_owner());
  EXPECT_FALSE(h.file_is_open());
  EXPECT_FALSE(h.group_is_open());
  EXPECT_EQ(h.file_id(), HDF5Handle::hid_default);
  EXPECT_EQ(h.group_id(), HDF5Handle::hid_default);
  EXPECT_NO_THROW(h.close_file());
  EXPECT_NO_THROW(h.close_group());
  EXPECT_EQ(h.open_file(HDF5Handle::Access::read_only), HDF5Handle::hid_default);
  EXPECT_EQ(h.open_file(HDF5Handle::Access::read_write), HDF5Handle::hid_default);
  EXPECT_EQ(h.open_group(), HDF5Handle::hid_default);
  EXPECT_EQ(h.file_name(), std::string{});
  EXPECT_EQ(h.group_name(), std::string{});
}

TEST(HDF5HandleDefault, open_file) {
  auto h = HDF5Handle();
  auto id = h.open_file(name_single_dataset, HDF5Handle::Access::read_only);
  ASSERT_NE(id, HDF5Handle::hid_default);
  ASSERT_FALSE(h.empty());
  EXPECT_EQ(h.file_name(), name_single_dataset);
  EXPECT_TRUE(h.file_owner());
  EXPECT_EQ(h.file_id(), id);
  EXPECT_FALSE(h.group_owner());
  EXPECT_EQ(h.group_id(), HDF5Handle::hid_default);
}

TEST(HDF5Handle, construct_from_file) {
  auto h = HDF5Handle(name_single_dataset);
  ASSERT_EQ(h.file_name(), name_single_dataset);
  EXPECT_TRUE(h.empty());
  EXPECT_FALSE(h.file_is_open());
  EXPECT_TRUE(h.file_owner());
}

TEST(HDF5Handle, construct_from_file_and_open_close) {
  auto h = HDF5Handle(name_single_dataset);
  EXPECT_TRUE(h.empty());
  EXPECT_TRUE(h.file_owner());
  auto id = h.open_file(HDF5Handle::Access::read_only);
  EXPECT_NE(id, HDF5Handle::hid_default);
  EXPECT_EQ(id, h.file_id());
  EXPECT_TRUE(h.file_owner());
  EXPECT_TRUE(h.file_is_open());
  EXPECT_NO_THROW(h.close_file());
  EXPECT_EQ(h.file_id(), HDF5Handle::hid_default);
  EXPECT_TRUE(h.file_owner());
}
