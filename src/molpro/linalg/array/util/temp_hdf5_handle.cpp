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
  bool set_erase_file_on_destroy_succeeded = handle.set_erase_file_on_destroy(true);
  assert(set_erase_file_on_destroy_succeeded);
  return handle;
}

std::string temp_group_name(const HDF5Handle &handle, const std::string &base_name) {
  if (handle.file_name().empty())
    throw std::runtime_error("temp_hdf5_handle_group: handle must have a file assigned");
  auto fid = handle.file_id();
  auto temp_handle = std::shared_ptr<HDF5Handle>{nullptr};
  if (not handle.file_is_open()) {
    temp_handle = std::make_shared<HDF5Handle>(handle.file_name());
    temp_handle->open_file(HDF5Handle::Access::read_only);
    fid = temp_handle->file_id();
  }
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
  auto res = hdf5_link_exists(fid, group_name);
  while (res > 0) {
    group_name = path + "_" + std::to_string(i++);
    res = hdf5_link_exists(fid, group_name);
  }
  if (res == -1)
    throw std::runtime_error("temp_hdf5_handle_group: failed to find a valid group name");
  return group_name;
}

HDF5Handle temp_hdf5_handle_group(const HDF5Handle &handle, const std::string &base_name) {
  if (not handle.file_owner())
    throw std::runtime_error(
        "temp_hdf5_handle_group: handle must be a file owner to create a temporary group in its copy");
  auto group_name = temp_group_name(handle, base_name);
  auto t = HDF5Handle(handle);
  t.open_file(HDF5Handle::Access::read_write);
  t.open_group(group_name);
  t.close_file();
  if (not t.set_erase_group_on_destroy(true))
    throw std::runtime_error("temp_hdf5_handle_group: failed to set_erase_group_on_destroy on the copy of handle");
  return t;
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro