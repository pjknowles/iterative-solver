#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data_util.h"

#include <molpro/linalg/array/HDF5Handle.h>

using molpro::linalg::array::util::file_exists;
using molpro::linalg::array::util::hdf5_link_exists;
using molpro::linalg::array::util::hdf5_open_file;
using molpro::linalg::array::util::HDF5Handle;

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

void check_dummy_state(const HDF5Handle& h) {
  EXPECT_TRUE(h.empty());
  EXPECT_FALSE(h.file_owner());
  EXPECT_FALSE(h.group_owner());
  EXPECT_FALSE(h.file_is_open());
  EXPECT_FALSE(h.group_is_open());
  EXPECT_EQ(h.file_id(), HDF5Handle::hid_default);
  EXPECT_EQ(h.group_id(), HDF5Handle::hid_default);
  EXPECT_TRUE(h.file_name().empty());
  EXPECT_TRUE(h.group_name().empty());
}

TEST_F(HDF5HandleDummyF, default_constructor) {
  check_dummy_state(*h);
  EXPECT_NO_THROW(h->close_file());
  EXPECT_NO_THROW(h->close_group());
  EXPECT_EQ(h->open_file(HDF5Handle::Access::read_only), HDF5Handle::hid_default);
  EXPECT_EQ(h->open_file(HDF5Handle::Access::read_write), HDF5Handle::hid_default);
  EXPECT_EQ(h->open_group(), HDF5Handle::hid_default);
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

TEST_F(HDF5HandleDummyF, open_file_with_name_check_access_type) {
  auto id = h->open_file(name_single_dataset, HDF5Handle::Access::read_only);
  ASSERT_TRUE(h->file_is_open());
  ASSERT_EQ(h->access_type(), HDF5Handle::Access::read_only);
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
  HDF5HandleOpenFileCreatF() : file_name{"test_new_file.hdf5"} { remove_file(); }
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

TEST_F(HDF5HandleOpenFileCreatF, create_file) {
  auto id = h->open_file(file_name, HDF5Handle::Access::read_write);
  EXPECT_TRUE(h->file_is_open());
  EXPECT_FALSE(h->empty());
  EXPECT_NE(id, HDF5Handle::hid_default);
}

TEST_F(HDF5HandleOpenFileCreatF, access_type) {
  EXPECT_EQ(h->access_type(), HDF5Handle::Access::none);
  auto id = h->open_file(file_name, HDF5Handle::Access::read_write);
  ASSERT_TRUE(h->file_is_open());
  EXPECT_EQ(h->access_type(), HDF5Handle::Access::read_write);
}

TEST_F(HDF5HandleOpenFileCreatF, create_group) {
  auto group_name = std::string{"/test-group1/test-group2"};
  h->open_file(file_name, HDF5Handle::Access::read_write);
  ASSERT_TRUE(h->file_is_open());
  auto gid = h->open_group(group_name);
  ASSERT_TRUE(h->group_is_open());
  ASSERT_EQ(h->group_id(), gid);
  ASSERT_NE(h->group_id(), HDF5Handle::hid_default);
  EXPECT_EQ(h->group_name(), group_name);
}

TEST_F(HDF5HandleDummyF, open_file_diff_name) {
  auto id_does_not_exist = h->open_file(name_inner_group_dataset + "-does-not-exist", HDF5Handle::Access::read_only);
  EXPECT_EQ(id_does_not_exist, HDF5Handle::hid_default)
      << "file does not exist and cannot be created with access type read_only";
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

TEST(HDF5Handle, construct_from_file_name) {
  auto h = HDF5Handle(name_single_dataset);
  ASSERT_EQ(h.file_name(), name_single_dataset);
  EXPECT_TRUE(h.empty());
  EXPECT_FALSE(h.file_is_open());
  EXPECT_TRUE(h.file_owner());
  auto fid = h.open_file(HDF5Handle::Access::read_only);
  ASSERT_TRUE(h.file_is_open());
  ASSERT_EQ(h.file_id(), fid);
  EXPECT_TRUE(h.file_owner());
  EXPECT_FALSE(h.empty());
  EXPECT_NE(fid, HDF5Handle::hid_default);
  h.close_file();
  ASSERT_FALSE(h.file_is_open());
  EXPECT_TRUE(h.file_owner());
  EXPECT_TRUE(h.empty());
  EXPECT_TRUE(H5Iis_valid(fid) == 0) << "file should have been closed";
}

TEST(HDF5Handle, construct_from_file_and_group_name) {
  std::string group_name = "/group1/group2";
  auto h = HDF5Handle(name_inner_group_dataset, group_name);
  ASSERT_EQ(h.file_name(), name_inner_group_dataset);
  ASSERT_EQ(h.group_name(), group_name);
  ASSERT_TRUE(h.empty());
  EXPECT_FALSE(h.file_is_open());
  EXPECT_TRUE(h.file_owner());
  EXPECT_TRUE(h.group_owner());
  auto gid = h.open_group();
  ASSERT_TRUE(h.group_is_open());
  ASSERT_TRUE(h.file_is_open());
  EXPECT_FALSE(h.empty());
  EXPECT_EQ(h.group_id(), gid);
  EXPECT_TRUE(H5Iis_valid(gid) > 0);
  EXPECT_TRUE(H5Iis_valid(h.file_id()) > 0);
  EXPECT_TRUE(h.file_owner());
  EXPECT_TRUE(h.group_owner());
  h.close_file();
  EXPECT_FALSE(h.file_is_open());
  EXPECT_FALSE(h.group_is_open()) << "closing the file should have closed the group as well";
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

class HDF5HandleUsingExistingObjectF : public ::testing::Test {
public:
  void SetUp() override {
    fid = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    gid = H5Gopen(fid, group_name.c_str(), H5P_DEFAULT);
    EXPECT_GT(H5Iis_valid(gid), 0);
    ASSERT_GT(H5Iis_valid(fid), 0);
  }
  void TearDown() override {
    if (!transferred_ownership_group)
      H5Gclose(gid);
    if (!transferred_ownership_file)
      H5Fclose(fid);
    EXPECT_EQ(H5Iis_valid(gid), 0);
    ASSERT_EQ(H5Iis_valid(fid), 0);
  }
  hid_t fid = -1;
  hid_t gid = -1;
  bool transferred_ownership_file = false;
  bool transferred_ownership_group = false;
  const std::string file_name = name_inner_group_dataset;
  const std::string group_name = "/group1";
  const std::string group_name2 = "/group1/group2";
};

TEST_F(HDF5HandleUsingExistingObjectF, construct_from_file_hid) {
  auto handle = HDF5Handle(fid);
  EXPECT_EQ(handle.file_name(), file_name);
  EXPECT_EQ(handle.file_id(), fid);
  EXPECT_TRUE(handle.file_is_open());
  EXPECT_FALSE(handle.file_owner());
  EXPECT_FALSE(handle.group_owner());
  EXPECT_FALSE(handle.group_is_open());
  EXPECT_TRUE(handle.group_name().empty());
  EXPECT_EQ(handle.group_id(), HDF5Handle::hid_default);
}

TEST_F(HDF5HandleUsingExistingObjectF, construct_from_file_hid_transfer_ownership) {
  transferred_ownership_file = true;
  auto handle = HDF5Handle(fid, transferred_ownership_file);
  EXPECT_EQ(handle.file_id(), fid);
  EXPECT_EQ(handle.file_name(), file_name);
  EXPECT_TRUE(handle.file_is_open());
  EXPECT_TRUE(handle.file_owner());
  EXPECT_FALSE(handle.group_owner());
  EXPECT_FALSE(handle.group_is_open());
  EXPECT_TRUE(handle.group_name().empty());
  EXPECT_EQ(handle.group_id(), HDF5Handle::hid_default);
}

TEST_F(HDF5HandleUsingExistingObjectF, construct_from_group_hid) {
  auto handle = HDF5Handle(gid);
  EXPECT_EQ(handle.group_id(), gid);
  EXPECT_EQ(handle.file_name(), file_name);
  EXPECT_TRUE(handle.group_name().find(group_name) != std::string::npos) << "handle stores absolute path";
  EXPECT_TRUE(handle.group_is_open());
  EXPECT_FALSE(handle.file_owner());
  EXPECT_FALSE(handle.group_owner());
  EXPECT_FALSE(handle.file_is_open());
  EXPECT_EQ(handle.file_id(), HDF5Handle::hid_default);
}

TEST_F(HDF5HandleUsingExistingObjectF, construct_from_group_hid_transfer_ownership) {
  transferred_ownership_group = true;
  auto handle = HDF5Handle(gid, transferred_ownership_group);
  EXPECT_EQ(handle.group_id(), gid);
  EXPECT_EQ(handle.file_name(), file_name);
  EXPECT_TRUE(handle.group_name().find(group_name) != std::string::npos) << "handle stores absolute path";
  EXPECT_TRUE(handle.group_is_open());
  EXPECT_TRUE(handle.file_owner());
  EXPECT_TRUE(handle.group_owner());
  EXPECT_FALSE(handle.file_is_open());
  EXPECT_EQ(handle.file_id(), HDF5Handle::hid_default);
}

TEST_F(HDF5HandleUsingExistingObjectF, open_and_close_file) {
  auto handle = HDF5Handle(fid);
  ASSERT_TRUE(handle.file_is_open());
  ASSERT_FALSE(handle.file_owner());
  auto id = handle.open_file(HDF5Handle::Access::read_only);
  ASSERT_TRUE(handle.file_is_open());
  ASSERT_FALSE(handle.file_owner());
  EXPECT_EQ(id, fid) << "file is already open, should return the same id back";
  handle.close_file();
  ASSERT_FALSE(handle.file_is_open());
  ASSERT_EQ(handle.file_id(), HDF5Handle::hid_default) << "without ownership close_file() resets to default id";
  EXPECT_GT(H5Iis_valid(fid), 0) << "without ownership close_file() does nothing to the original object";
}

TEST_F(HDF5HandleUsingExistingObjectF, open_file_wrong_access_type) {
  auto handle = HDF5Handle(fid);
  ASSERT_TRUE(handle.file_is_open());
  ASSERT_FALSE(handle.file_owner());
  auto id = handle.open_file(HDF5Handle::Access::read_write);
  ASSERT_TRUE(handle.file_is_open());
  ASSERT_FALSE(handle.file_owner());
  EXPECT_EQ(id, HDF5Handle::hid_default)
      << "attempt to open with different access type, returns default id to indicate an error";
  EXPECT_EQ(handle.file_id(), fid) << "stored id should not be modified";
}

TEST_F(HDF5HandleUsingExistingObjectF, open_file_no_ownership_and_was_closed_outside) {
  auto handle = HDF5Handle(fid);
  ASSERT_TRUE(handle.file_is_open());
  ASSERT_FALSE(handle.file_owner());
  transferred_ownership_file = true;
  H5Fclose(fid);
  ASSERT_FALSE(handle.file_is_open());
  auto id = handle.open_file(HDF5Handle::Access::read_only);
  ASSERT_FALSE(handle.file_is_open());
  ASSERT_FALSE(handle.file_owner());
  EXPECT_EQ(id, HDF5Handle::hid_default) << "handle is not owning and file was closed outside";
  EXPECT_EQ(handle.file_id(), fid) << "stored id should not be modified";
  id = handle.open_file(file_name, HDF5Handle::Access::read_only);
  ASSERT_FALSE(handle.file_is_open());
  ASSERT_FALSE(handle.file_owner());
  EXPECT_EQ(id, HDF5Handle::hid_default) << "handle is not owning so cannot reopen the file";
}

TEST_F(HDF5HandleUsingExistingObjectF, open_and_close_file_with_ownership) {
  transferred_ownership_file = true;
  auto handle = HDF5Handle(fid, transferred_ownership_file);
  ASSERT_TRUE(handle.file_is_open());
  ASSERT_TRUE(handle.file_owner());
  ASSERT_EQ(handle.file_id(), fid);
  auto id = handle.open_file(HDF5Handle::Access::read_only);
  ASSERT_TRUE(handle.file_is_open());
  ASSERT_TRUE(handle.file_owner());
  ASSERT_EQ(id, fid) << "file was already open";
  handle.close_file();
  ASSERT_FALSE(handle.file_is_open()) << "with ownership should be able to close the file";
  ASSERT_TRUE(handle.file_owner());
}

TEST_F(HDF5HandleUsingExistingObjectF, construct_from_file_no_own__open_and_close_group) {
  auto handle = HDF5Handle(fid);
  ASSERT_TRUE(handle.file_is_open());
  ASSERT_FALSE(handle.file_owner());
  auto gid = handle.open_group();
  ASSERT_FALSE(handle.group_is_open()) << "no group was assigned";
  ASSERT_EQ(gid, HDF5Handle::hid_default);
  gid = handle.open_group(group_name);
  ASSERT_TRUE(handle.group_is_open()) << "can open the group even without owning the file";
  ASSERT_NE(gid, HDF5Handle::hid_default);
  handle.close_group();
  ASSERT_FALSE(handle.group_is_open()) << "opening group gave it ownership of the group, can close the file";
}

TEST_F(HDF5HandleUsingExistingObjectF, construct_from_group_no_own__open_and_close_group) {
  auto handle = HDF5Handle(gid);
  ASSERT_FALSE(handle.file_is_open());
  ASSERT_FALSE(handle.file_owner());
  ASSERT_TRUE(handle.group_is_open());
  ASSERT_FALSE(handle.group_owner());
  auto id = handle.open_group();
  EXPECT_EQ(id, gid) << "group is already open, just return the id";
  handle.close_group();
  ASSERT_FALSE(handle.group_is_open());
  ASSERT_EQ(handle.group_id(), HDF5Handle::hid_default) << "without ownership close_group() resets to default id";
  EXPECT_GT(H5Iis_valid(gid), 0) << "without ownership close_group() does nothing to the original object";
  id = handle.open_group(group_name2);
  EXPECT_TRUE(handle.group_is_open());
  EXPECT_TRUE(handle.group_owner());
  EXPECT_NE(id, HDF5Handle::hid_default) << "can reassign the group even without ownership";
  EXPECT_NE(id, gid);
  EXPECT_GT(H5Iis_valid(gid), 0) << "without ownership reassignment does not close the original group";
  EXPECT_TRUE(handle.file_is_open()) << "opens the file";
  EXPECT_TRUE(handle.file_owner()) << "gives ownership to the file";
  unsigned intent = 0;
  auto err = H5Fget_intent(handle.file_id(), &intent);
  EXPECT_EQ(intent, H5F_ACC_RDONLY) << "the file should be read-only";
  handle.close_group();
  EXPECT_FALSE(handle.group_is_open());
  EXPECT_GT(H5Iis_valid(gid), 0) << "without ownership reassignment does not close the original group";
}

TEST_F(HDF5HandleUsingExistingObjectF, construct_from_group_no_own__open_file) {
  auto handle = HDF5Handle(gid);
  ASSERT_FALSE(handle.file_owner());
  auto id = handle.open_file(HDF5Handle::Access::read_only);
  ASSERT_FALSE(handle.file_is_open()) << "without ownership should not be able to open the file";
  EXPECT_EQ(id, HDF5Handle::hid_default);
  id = handle.open_file(file_name + "-does-not-exist-2", HDF5Handle::Access::read_write);
  ASSERT_FALSE(handle.file_is_open()) << "a different file is already assigned";
  EXPECT_EQ(id, HDF5Handle::hid_default);
}

TEST_F(HDF5HandleUsingExistingObjectF, construct_from_group_with_own__open_file) {
  transferred_ownership_group = true;
  auto handle = HDF5Handle(gid, transferred_ownership_group);
  ASSERT_TRUE(handle.file_owner());
  ASSERT_FALSE(handle.file_is_open());
  auto id = handle.open_file(HDF5Handle::Access::read_only);
  ASSERT_TRUE(handle.file_is_open()) << "with ownership can open the underlying file";
  EXPECT_NE(id, HDF5Handle::hid_default);
  handle.close_file();
  ASSERT_FALSE(handle.file_is_open());
  id = handle.open_file(file_name + "-does-not-exist-2", HDF5Handle::Access::read_write);
  ASSERT_FALSE(handle.file_is_open()) << "reassignment of files is forbidden";
  EXPECT_EQ(id, HDF5Handle::hid_default);
}

inline void copy_constructor_same_state(const HDF5Handle& h, const HDF5Handle& hsource) {
  EXPECT_EQ(h.access_type(), hsource.access_type());
  EXPECT_EQ(h.file_name(), hsource.file_name());
  EXPECT_EQ(h.group_name(), hsource.group_name());
  EXPECT_EQ(h.file_is_open(), hsource.file_is_open());
  EXPECT_EQ(h.group_is_open(), hsource.group_is_open());
  EXPECT_EQ(h.file_owner(), hsource.file_owner());
  EXPECT_EQ(h.group_owner(), hsource.group_owner());
  EXPECT_EQ(h.empty(), hsource.empty());
}

TEST(HDF5Handle, copy_constructor_source_from_dummy) {
  auto hsource = HDF5Handle();
  auto h = HDF5Handle(hsource);
  check_dummy_state(h);
  copy_constructor_same_state(h, hsource);
  EXPECT_EQ(h.file_id(), hsource.file_id());
  EXPECT_EQ(h.group_id(), hsource.group_id());
}

TEST(HDF5Handle, copy_constructor_source_from_file_name_file_closed) {
  auto hsource = HDF5Handle(name_inner_group_dataset);
  auto h = HDF5Handle(hsource);
  copy_constructor_same_state(h, hsource);
  EXPECT_EQ(h.file_id(), hsource.file_id());
  EXPECT_EQ(h.group_id(), hsource.group_id());
}

TEST(HDF5Handle, copy_constructor_source_from_file_name_file_open) {
  auto hsource = HDF5Handle(name_inner_group_dataset);
  hsource.open_file(HDF5Handle::Access::read_only);
  hsource.open_group("group1");
  ASSERT_TRUE(hsource.file_is_open());
  ASSERT_TRUE(hsource.group_is_open());
  auto h = HDF5Handle(hsource);
  copy_constructor_same_state(h, hsource);
  EXPECT_NE(h.file_id(), hsource.file_id());
  EXPECT_NE(h.group_id(), hsource.group_id());
}

TEST(HDF5Handle, copy_constructor_source_from_file_and_group_name) {
  auto hsource = HDF5Handle(name_inner_group_dataset, "group1");
  auto h = HDF5Handle(hsource);
  copy_constructor_same_state(h, hsource);
  EXPECT_EQ(h.file_id(), hsource.file_id());
  EXPECT_EQ(h.group_id(), hsource.group_id());
}

TEST(HDF5Handle, copy_assignment_source_from_file_and_group_name_to_diff_file_name) {
  auto hsource = HDF5Handle(name_inner_group_dataset, "group1");
  auto h = HDF5Handle(name_single_dataset);
  hsource = h;
  copy_constructor_same_state(h, hsource);
  EXPECT_EQ(h.file_id(), hsource.file_id());
  EXPECT_EQ(h.group_id(), hsource.group_id());
}

TEST_F(HDF5HandleUsingExistingObjectF, copy_constructor_source_from_file_obj) {
  auto hsource = HDF5Handle(fid);
  auto h = HDF5Handle(hsource);
  copy_constructor_same_state(h, hsource);
  EXPECT_EQ(h.file_id(), hsource.file_id());
  EXPECT_EQ(h.group_id(), hsource.group_id());
}

TEST_F(HDF5HandleUsingExistingObjectF, copy_constructor_source_from_file_obj_with_own) {
  transferred_ownership_file = true;
  auto hsource = HDF5Handle(fid, transferred_ownership_file);
  auto h = HDF5Handle(hsource);
  copy_constructor_same_state(h, hsource);
  EXPECT_NE(h.file_id(), hsource.file_id());
  EXPECT_EQ(h.group_id(), hsource.group_id());
}

TEST_F(HDF5HandleUsingExistingObjectF, copy_constructor_source_from_group_obj) {
  auto hsource = HDF5Handle(gid);
  auto h = HDF5Handle(hsource);
  copy_constructor_same_state(h, hsource);
  EXPECT_EQ(h.file_id(), hsource.file_id());
  EXPECT_EQ(h.group_id(), hsource.group_id());
}

inline void copy_constructor_same_state_grom_group_obj(const HDF5Handle& h, const HDF5Handle& hsource) {
  EXPECT_EQ(h.access_type(), HDF5Handle::Access::read_only)
      << "even though source had a closed file, opening a group automatically opens the file";
  EXPECT_EQ(h.file_name(), hsource.file_name());
  EXPECT_EQ(h.group_name(), hsource.group_name());
  EXPECT_TRUE(h.file_is_open()) << "even though source had a closed file, opening a group automatically opens the file";
  EXPECT_EQ(h.group_is_open(), hsource.group_is_open());
  EXPECT_EQ(h.file_owner(), hsource.file_owner());
  EXPECT_EQ(h.group_owner(), hsource.group_owner());
  EXPECT_EQ(h.empty(), hsource.empty());
  EXPECT_NE(h.file_id(), hsource.file_id());
  EXPECT_NE(h.group_id(), hsource.group_id());
}

TEST_F(HDF5HandleUsingExistingObjectF, copy_constructor_source_from_group_obj_with_own) {
  transferred_ownership_group = true;
  auto hsource = HDF5Handle(gid, transferred_ownership_group);
  auto h = HDF5Handle(hsource);
  copy_constructor_same_state_grom_group_obj(h, hsource);
}

TEST_F(HDF5HandleUsingExistingObjectF, assignment_operator_from_file_name_open_to_group_obj_with_own) {
  transferred_ownership_group = true;
  auto hsource = HDF5Handle(gid, transferred_ownership_group);
  ASSERT_NE(file_name, name_single_dataset);
  auto h = HDF5Handle(name_single_dataset);
  auto id = h.open_file(HDF5Handle::Access::read_only);
  h = hsource;
  copy_constructor_same_state_grom_group_obj(h, hsource);
  EXPECT_NE(h.file_id(), id);
  EXPECT_NE(h.file_name(), name_single_dataset);
}

TEST(HDF5Handle, move_assignment_source_from_dummy) {
  HDF5Handle&& hsource = HDF5Handle();
  auto h = HDF5Handle();
  h = std::forward<HDF5Handle>(hsource);
  check_dummy_state(h);
  check_dummy_state(hsource);
}

TEST(HDF5Handle, move_assignment_source_from_file_group_name_open) {
  HDF5Handle&& hsource = HDF5Handle(name_single_dataset, "/");
  hsource.open_group();
  ASSERT_TRUE(hsource.file_is_open());
  ASSERT_TRUE(hsource.group_is_open());
  auto fid_source = hsource.file_id();
  auto gid_source = hsource.group_id();
  auto h = HDF5Handle(name_inner_group_dataset, "group1");
  h.open_group();
  ASSERT_TRUE(h.file_is_open());
  ASSERT_TRUE(h.group_is_open());
  auto fid = h.file_id();
  auto gid = h.group_id();
  h = std::forward<HDF5Handle>(hsource);
  check_dummy_state(hsource);
  EXPECT_TRUE(H5Iis_valid(fid) == 0) << "previous file should have been closed";
  EXPECT_TRUE(H5Iis_valid(gid) == 0) << "previous group should have been closed";
  EXPECT_EQ(h.file_id(), fid_source) << "take ownership of the previous id";
  EXPECT_EQ(h.group_id(), gid_source) << "take ownership of the previous id";
  EXPECT_TRUE(h.file_is_open());
  EXPECT_TRUE(h.group_is_open());
  EXPECT_TRUE(h.file_owner());
  EXPECT_TRUE(h.group_owner());
  EXPECT_EQ(h.file_name(), name_single_dataset);
  EXPECT_EQ(h.group_name(), "/");
}

TEST(HDF5Handle, move_constructor_from_file_group_name_open) {
  auto file_name = name_inner_group_dataset;
  auto group_name = "/group1";
  HDF5Handle&& hsource = HDF5Handle(file_name, group_name);
  hsource.open_group();
  ASSERT_TRUE(hsource.file_is_open());
  ASSERT_TRUE(hsource.group_is_open());
  auto fid_source = hsource.file_id();
  auto gid_source = hsource.group_id();
  auto h = HDF5Handle(std::forward<HDF5Handle>(hsource));
  check_dummy_state(hsource);
  EXPECT_TRUE(h.file_is_open());
  EXPECT_TRUE(h.group_is_open());
  EXPECT_EQ(h.file_id(), fid_source) << "take ownership of the previous id";
  EXPECT_EQ(h.group_id(), gid_source) << "take ownership of the previous id";
  EXPECT_TRUE(h.file_owner());
  EXPECT_TRUE(h.group_owner());
  EXPECT_EQ(h.file_name(), file_name);
  EXPECT_EQ(h.group_name(), group_name);
}

TEST_F(HDF5HandleUsingExistingObjectF, move_assignment_operator_from_group_obj_with_own_to_file_name) {
  transferred_ownership_group = true;
  HDF5Handle&& hsource = HDF5Handle(gid, transferred_ownership_group);
  hsource.open_group();
  ASSERT_FALSE(hsource.file_is_open());
  ASSERT_TRUE(hsource.group_is_open());
  auto fid_source = hsource.file_id();
  auto gid_source = hsource.group_id();
  auto h = HDF5Handle(name_single_dataset);
  h.open_file(HDF5Handle::Access::read_only);
  ASSERT_TRUE(h.file_is_open());
  ASSERT_FALSE(h.group_is_open());
  auto fid_target = h.file_id();
  h = std::forward<HDF5Handle>(hsource);
  check_dummy_state(hsource);
  EXPECT_TRUE(H5Iis_valid(fid_target) == 0) << "previous file should have been closed";
  EXPECT_EQ(h.file_id(), fid_source) << "take ownership of the previous id";
  EXPECT_EQ(h.group_id(), gid_source) << "take ownership of the previous id";
  EXPECT_FALSE(h.file_is_open());
  EXPECT_TRUE(h.group_is_open());
  EXPECT_TRUE(h.file_owner());
  EXPECT_TRUE(h.group_owner());
  EXPECT_EQ(h.file_name(), file_name);
  EXPECT_EQ(h.group_name(), group_name);
}

TEST_F(HDF5HandleUsingExistingObjectF, move_constructor_operator_from_group_obj_with_own) {
  transferred_ownership_group = true;
  HDF5Handle&& hsource = HDF5Handle(gid, transferred_ownership_group);
  hsource.open_group();
  ASSERT_FALSE(hsource.file_is_open());
  ASSERT_TRUE(hsource.group_is_open());
  auto fid_source = hsource.file_id();
  auto gid_source = hsource.group_id();
  auto h = HDF5Handle(std::forward<HDF5Handle>(hsource));
  check_dummy_state(hsource);
  EXPECT_FALSE(h.file_is_open());
  EXPECT_TRUE(h.group_is_open());
  EXPECT_EQ(h.file_id(), fid_source) << "take ownership of the previous id";
  EXPECT_EQ(h.group_id(), gid_source) << "take ownership of the previous id";
  EXPECT_TRUE(h.file_owner());
  EXPECT_TRUE(h.group_owner());
  EXPECT_EQ(h.file_name(), file_name);
  EXPECT_EQ(h.group_name(), group_name);
}
