#include "HDF5Handle.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

const hid_t HDF5Handle::hid_default;

HDF5Handle::HDF5Handle(std::string file) : m_file_name{std::move(file)}, m_file_owner{true} {}
HDF5Handle::HDF5Handle(std::string file, std::string group)
    : m_file_name{std::move(file)}, m_group_name{std::move(group)}, m_file_owner{true}, m_group_owner{true} {}
HDF5Handle::HDF5Handle(hid_t hid, bool transfer_ownership) {
  if (H5Iis_valid(hid) > 0) {
    auto type = H5Iget_type(hid);
    m_file_name = hdf5_get_file_name(hid);
    m_file_owner = transfer_ownership;
    if (type == H5I_FILE) {
      m_file_hid = hid;
    } else if (type == H5I_GROUP) {
      m_group_hid = hid;
      m_group_owner = transfer_ownership;
      m_group_name = hdf5_get_object_name(m_group_hid);
    }
  }
}

HDF5Handle::HDF5Handle(const HDF5Handle &source) { *this = source; }

HDF5Handle &HDF5Handle::operator=(const HDF5Handle &source) {
  if (file_is_open())
    close_file();
  m_file_hid = source.m_file_hid;
  m_group_hid = source.m_group_hid;
  m_file_name = source.m_file_name;
  m_group_name = source.m_group_name;
  m_file_owner = source.m_file_owner;
  m_group_owner = source.m_group_owner;
  m_erase_file_on_destroy = source.m_erase_file_on_destroy;
  if (m_file_owner) {
    m_file_hid = hid_default;
    if (source.file_is_open())
      HDF5Handle::open_file(source.access_type());
  }
  if (m_group_owner) {
    m_group_hid = hid_default;
    if (source.group_is_open())
      HDF5Handle::open_group();
  }
  return *this;
}

HDF5Handle::HDF5Handle(HDF5Handle &&source) noexcept : HDF5Handle{} { *this = std::forward<HDF5Handle>(source); }

HDF5Handle &HDF5Handle::operator=(HDF5Handle &&source) noexcept {
  if (file_is_open())
    close_file();
  m_file_hid = source.m_file_hid;
  m_group_hid = source.m_group_hid;
  m_file_name = source.m_file_name;
  m_group_name = source.m_group_name;
  m_file_owner = source.m_file_owner;
  m_group_owner = source.m_group_owner;
  m_erase_file_on_destroy = source.m_erase_file_on_destroy;
  source.m_file_owner = false;
  source.m_group_owner = false;
  source.m_erase_file_on_destroy = false;
  auto dummy = HDF5Handle{};
  source = dummy;
  return *this;
}

HDF5Handle::~HDF5Handle() {
  HDF5Handle::close_file();
  if (m_erase_file_on_destroy) {
    if (file_exists(file_name())) {
      std::remove(file_name().c_str());
    }
  } else if (m_erase_group_on_destroy) {
    HDF5Handle::open_file(Access::read_write);
    HDF5Handle::open_group();
    H5Ldelete(file_id(), group_name().c_str(), H5P_DEFAULT);
  }
}
hid_t HDF5Handle::open_file(HDF5Handle::Access type) {
  if (file_is_open()) {
    if (access_type() == type)
      return m_file_hid;
    else
      return hid_default;
  }
  if (!m_file_owner || m_file_name.empty())
    return hid_default;
  auto plist = _open_plist();
  if (file_exists(m_file_name)) {
    if (type == Access::read_only)
      m_file_hid = H5Fopen(m_file_name.c_str(), H5F_ACC_RDONLY, plist);
    else
      m_file_hid = H5Fopen(m_file_name.c_str(), H5F_ACC_RDWR, plist);
  } else {
    if (type == Access::read_only)
      m_file_hid = hid_default;
    else
      m_file_hid = H5Fcreate(m_file_name.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, plist);
  }
  if (plist != H5P_DEFAULT)
    H5Pclose(plist);
  if (m_file_hid < 0) {
    m_file_hid = hid_default;
  }
  return m_file_hid;
}
hid_t HDF5Handle::open_file(const std::string &file, HDF5Handle::Access type) {
  if (!m_file_name.empty()) { // file has been assigned
    if (m_file_name != file)  // attempting to reassign
      return hid_default;
  } else {
    m_file_name = file;  // first assignment
    m_file_owner = true; // takes ownership
  }
  return open_file(type);
}
void HDF5Handle::close_file() {
  close_group();
  if (!m_file_owner) {
    m_file_hid = hid_default;
    return;
  }
  if (file_is_open())
    H5Fclose(m_file_hid);
  m_file_hid = hid_default;
}
bool HDF5Handle::file_is_open() const {
  if (m_file_hid == hid_default)
    return false;
  return H5Iis_valid(m_file_hid) > 0;
}
bool HDF5Handle::group_is_open() const {
  if (m_group_hid == hid_default)
    return false;
  return H5Iis_valid(m_group_hid) > 0;
}
std::string HDF5Handle::file_name() const { return m_file_name; }
std::string HDF5Handle::group_name() const { return m_group_name; }
void HDF5Handle::close_group() {
  if (!m_group_owner) {
    m_group_hid = hid_default;
    return;
  }
  if (group_is_open())
    H5Gclose(m_group_hid);
  m_group_hid = hid_default;
}
hid_t HDF5Handle::open_group() {
  if (group_is_open())
    return m_group_hid;
  if (!m_group_owner || m_group_name.empty())
    return hid_default;
  if (!file_is_open()) {
    auto prev_ownership = m_file_owner;
    m_file_owner = true;
    auto fid = open_file(Access::read_only);
    if (fid == hid_default) { // failed to open the file, give back previous ownership
      m_file_owner = prev_ownership;
      return hid_default;
    }
  }

  if (hdf5_link_exists(m_file_hid, m_group_name) < 1) {
    auto lcpl = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(lcpl, 1);
    m_group_hid = H5Gcreate(m_file_hid, m_group_name.c_str(), lcpl, H5P_DEFAULT, H5P_DEFAULT);
    H5Pclose(lcpl);
  } else {
    m_group_hid = H5Gopen(m_file_hid, m_group_name.c_str(), H5P_DEFAULT);
  }
  if (m_group_hid < 0) {
    m_group_hid = hid_default;
  }
  return m_group_hid;
}
hid_t HDF5Handle::open_group(const std::string &group) {
  if (m_group_name != group)
    close_group();
  m_group_name = group;
  m_group_owner = true;
  return open_group();
}
HDF5Handle::Access HDF5Handle::access_type() const {
  if (!file_is_open())
    return Access::none;
  unsigned int intent;
  auto err = H5Fget_intent(m_file_hid, &intent);
  if (err < 0)
    return Access::unknown;
  if (H5F_ACC_RDONLY == intent)
    return Access::read_only;
  else if (H5F_ACC_RDWR & intent)
    return Access::read_write;
  return Access::unknown;
}
hid_t HDF5Handle::file_id() const { return m_file_hid; }
hid_t HDF5Handle::group_id() const { return m_group_hid; }
bool HDF5Handle::file_owner() const { return m_file_owner; }
bool HDF5Handle::group_owner() const { return m_group_owner; }
bool HDF5Handle::empty() const { return m_file_hid == hid_default && m_group_hid == hid_default; }
hid_t HDF5Handle::_open_plist() { return H5P_DEFAULT; }

bool HDF5Handle::set_erase_file_on_destroy(bool value) {
  if (value && !erasable())
    return false;
  m_erase_file_on_destroy = value;
  return true;
}

bool HDF5Handle::set_erase_group_on_destroy(bool value) {
  if (value && !m_group_owner)
    return false;
  m_erase_group_on_destroy = value;
  return true;
}

bool HDF5Handle::erasable() {
  bool file_owner = m_file_owner;
  bool group_owner = (m_group_hid == hid_default) || m_group_owner;
  return file_owner && group_owner;
}

bool file_exists(const std::string &fname) { return !std::ifstream{fname}.fail(); }

hid_t hdf5_open_file(const std::string &fname, bool read_only) {
  hid_t id;
  if (read_only)
    id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  else
    id = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  return id;
}

bool hdf5_file_is_open(hid_t file_id) { return H5Fget_obj_count(file_id, H5F_OBJ_FILE) > 0; }

std::string hdf5_get_object_name(hid_t id) {
  auto size = H5Iget_name(id, nullptr, 0);
  std::string result;
  result.resize(size + 1);
  H5Iget_name(id, &result[0], size + 1);
  result.resize(size);
  return result;
}
std::string hdf5_get_file_name(hid_t id) {
  auto size = H5Fget_name(id, nullptr, 0);
  std::string result;
  result.resize(size + 1);
  H5Fget_name(id, &result[0], size + 1);
  result.resize(size); // get rid of the null terminator
  return result;
}
//  "/group1/group2/dataset/does_not_exist"
htri_t hdf5_link_exists(hid_t id, std::string path) {
  path.erase(path.begin(), std::find_if(path.begin(), path.end(), [](auto &el) { return !std::isspace(el); }));
  path.erase(std::find_if(path.rbegin(), path.rend(), [](auto &el) { return !std::isspace(el); }).base(), path.end());
  if (path.empty())
    return -1;
  htri_t res = -1;
  size_t n = 0;
  while (n != std::string::npos) {
    n = path.find('/', n);
    if (n != std::string::npos)
      ++n;
    auto sub = path.substr(0, n);
    res = H5Lexists(id, sub.c_str(), H5P_DEFAULT);
    if (res < 1)
      break;
  }
  return res;
}

template struct TempHandle<HDF5Handle>;
} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro
