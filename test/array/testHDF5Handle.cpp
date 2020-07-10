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

class HDF5HandleDummyF : public ::testing::Test {
public:
  HDF5HandleDummyF() = default;
  void SetUp() override { h = std::make_unique<HDF5Handle>(); }
  void TearDown() override { h.reset(); }
  std::unique_ptr<HDF5Handle> h;
};

TEST_F(HDF5HandleDummyF, default_constructor) {
  ASSERT_TRUE(h->empty());
  EXPECT_FALSE(h->file_owner());
  EXPECT_FALSE(h->group_owner());
  EXPECT_FALSE(h->file_is_open());
  EXPECT_FALSE(h->group_is_open());
  EXPECT_EQ(h->file_id(), HDF5Handle::hid_default);
  EXPECT_EQ(h->group_id(), HDF5Handle::hid_default);
  EXPECT_NO_THROW(h->close_file());
  EXPECT_NO_THROW(h->close_group());
  EXPECT_EQ(h->open_file(HDF5Handle::Access::read_only), HDF5Handle::hid_default);
  EXPECT_EQ(h->open_file(HDF5Handle::Access::read_write), HDF5Handle::hid_default);
  EXPECT_EQ(h->open_group(), HDF5Handle::hid_default);
  EXPECT_EQ(h->file_name(), std::string{});
  EXPECT_EQ(h->group_name(), std::string{});
}

TEST_F(HDF5HandleDummyF, open_file) {
  auto id = h->open_file(HDF5Handle::Access::read_only);
  ASSERT_EQ(id, HDF5Handle::hid_default) << "dummy handle has no file assigned to it";
  EXPECT_TRUE(h->empty());
}

TEST_F(HDF5HandleDummyF, open_file_with_name) {
  auto id = h->open_file(name_single_dataset, HDF5Handle::Access::read_only);
  ASSERT_NE(id, HDF5Handle::hid_default);
  ASSERT_FALSE(h->empty());
  EXPECT_EQ(h->file_name(), name_single_dataset);
  EXPECT_TRUE(h->file_owner());
  EXPECT_EQ(h->file_id(), id);
  EXPECT_FALSE(h->group_owner());
  EXPECT_EQ(h->group_id(), HDF5Handle::hid_default);
}

TEST_F(HDF5HandleDummyF, open_file_repeat) {
  auto id_first = h->open_file(name_single_dataset, HDF5Handle::Access::read_only);
  ASSERT_NE(id_first, HDF5Handle::hid_default);
  ASSERT_FALSE(h->empty());
  auto id2 = h->open_file(HDF5Handle::Access::read_only);
  EXPECT_EQ(id2, id_first) << "opening an open file with the same access should return previous id";
  auto id3 = h->open_file(HDF5Handle::Access::read_write);
  EXPECT_NE(id3, id_first) << "opening an opened file with different access should return hid_default";
  EXPECT_EQ(id3, HDF5Handle::hid_default) << "opening an opened file with different access should return hid_default";
  auto id_diff_file = h->open_file(name_inner_group_dataset, HDF5Handle::Access::read_only);
  EXPECT_EQ(id_diff_file, HDF5Handle::hid_default)
      << "a file with different name was already assigned, should return hid_default";
}

TEST_F(HDF5HandleDummyF, open_close_group) {
  h->open_file(name_inner_group_dataset, HDF5Handle::Access::read_only);
  ASSERT_TRUE(h->file_is_open());
  auto gid = h->open_group("group1");
  EXPECT_NE(gid, HDF5Handle::hid_default);
  ASSERT_TRUE(h->group_is_open());
  EXPECT_EQ(h->group_name(), "group1");
  EXPECT_TRUE(h->group_owner());
  h->close_group();
  ASSERT_FALSE(h->group_is_open());
  EXPECT_EQ(h->group_name(), "group1");
  EXPECT_TRUE(h->group_owner());
  h->open_group("group1");
  ASSERT_TRUE(h->group_is_open());
}

TEST_F(HDF5HandleDummyF, open_close_group_nested) {
  h->open_file(name_inner_group_dataset, HDF5Handle::Access::read_only);
  ASSERT_TRUE(h->file_is_open());
  auto gid = h->open_group("group1/group2");
  EXPECT_NE(gid, HDF5Handle::hid_default);
  ASSERT_TRUE(h->group_is_open());
  EXPECT_EQ(h->group_name(), "group1/group2");
  EXPECT_TRUE(h->group_owner());
  h->close_group();
  ASSERT_FALSE(h->group_is_open());
  EXPECT_EQ(h->group_name(), "group1/group2");
  EXPECT_TRUE(h->group_owner());
}

TEST_F(HDF5HandleDummyF, open_close_group_reassignment) {
  h->open_file(name_inner_group_dataset, HDF5Handle::Access::read_only);
  ASSERT_TRUE(h->file_is_open());
  auto gid_first = h->open_group("group1");
  EXPECT_NE(gid_first, HDF5Handle::hid_default);
  ASSERT_TRUE(h->group_is_open());
  ASSERT_TRUE(h->group_owner());
  auto gid = h->open_group("group1/group2");
  ASSERT_TRUE(h->group_is_open());
  EXPECT_NE(gid, HDF5Handle::hid_default);
  EXPECT_EQ(h->group_name(), "group1/group2");
  EXPECT_TRUE(h->group_owner());
  EXPECT_TRUE(H5Iis_valid(gid_first) == 0) << "previous group should have been closed";
}

struct HDF5HandleOpenFileCreatF : public HDF5HandleDummyF {
  HDF5HandleOpenFileCreatF() : file_name{name_inner_group_dataset + "-does-not-exist"} { remove_file(); }
  ~HDF5HandleOpenFileCreatF() { remove_file(); }
  void remove_file() {
    if (file_exists(file_name))
      std::remove(file_name.c_str());
  }
  const std::string file_name;
};

TEST_F(HDF5HandleOpenFileCreatF, wrong_access) {
  auto id = h->open_file(file_name, HDF5Handle::Access::read_only);
  EXPECT_EQ(id, HDF5Handle::hid_default) << "file does not exist and cannot be crated with access type read_only";
  EXPECT_FALSE(h->file_is_open());
  EXPECT_TRUE(h->empty());
}

TEST_F(HDF5HandleOpenFileCreatF, create) {
  auto id = h->open_file(file_name, HDF5Handle::Access::read_write);
  EXPECT_TRUE(h->file_is_open());
  EXPECT_FALSE(h->empty());
  EXPECT_NE(id, HDF5Handle::hid_default);
}

TEST_F(HDF5HandleOpenFileCreatF, open_group) {}

TEST_F(HDF5HandleDummyF, open_file_diff_name) {
  auto id_does_not_exist = h->open_file(name_inner_group_dataset + "-does-not-exist", HDF5Handle::Access::read_only);
  EXPECT_EQ(id_does_not_exist, HDF5Handle::hid_default)
      << "file does not exist and cannot be crated with access type read_only";
}

TEST_F(HDF5HandleDummyF, file_is_open) {
  ASSERT_FALSE(h->file_is_open());
  auto id = h->open_file(name_single_dataset, HDF5Handle::Access::read_only);
  ASSERT_TRUE(h->file_is_open());
  H5Fclose(id);
  ASSERT_FALSE(h->file_is_open()) << "handle should check if the file has been closed outside";
  EXPECT_NO_THROW(h->open_file(name_single_dataset, HDF5Handle::Access::read_only));
  ASSERT_NE(h->file_id(), HDF5Handle::hid_default) << "if file was closed outside, it can be opened by the handle";
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
