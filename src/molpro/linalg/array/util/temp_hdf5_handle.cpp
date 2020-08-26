#include "temp_hdf5_handle.h"
#include "temp_file.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <stdexcept>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

HDF5Handle temp_hdf5_handle(const std::string &base_name) {
  auto fname = temp_file_name(base_name, ".hdf5");
  auto handle = HDF5Handle(fname);
  assert(handle.set_erase_on_destroy(true));
  return handle;
}

std::string temp_group_name(const HDF5Handle &handle, const std::string &base_name) {
  if (handle.file_name().empty())
    throw std::runtime_error("temp_hdf5_handle_group: handle must have a file assigned");
  if (handle.file_id() == HDF5Handle::hid_default)
    throw std::runtime_error("temp_hdf5_handle_group: handle must have a file opened");
  auto path = base_name;
  path.erase(path.begin(), std::find_if(path.begin(), path.end(), [](auto &el) { return !std::isspace(el); }));
  path.erase(std::find_if(path.rbegin(), path.rend(), [](auto &el) { return !std::isspace(el); }).base(), path.end());
  if (path.empty())
    path = "/";
  bool is_absolute_path = (path.rfind('/', 0) == 0);
  if (not is_absolute_path) {
    if (handle.group_name().empty())
      throw std::runtime_error(
          "temp_hdf5_handle_group: if handle does not have a group_name than base_name must be an absolute path");
    path = handle.group_name() + "/" + path;
  }
  size_t i = 0;
  auto group_name = path;
  auto res = hdf5_link_exists(handle.file_id(), group_name);
  while (res > 0) {
    group_name = path + "_" + std::to_string(i++);
    res = hdf5_link_exists(handle.file_id(), group_name);
  }
  if (res == -1)
    throw std::runtime_error("temp_hdf5_handle_group: failed to find a valid group name");
  return group_name;
}

HDF5Handle temp_hdf5_handle_group(const HDF5Handle &handle, const std::string &base_name) {
  auto t = HDF5Handle(handle);
  auto group_name = temp_group_name(handle, base_name);
  t.open_group(group_name);
  return t;
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro