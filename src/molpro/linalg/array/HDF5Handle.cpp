#include "HDF5Handle.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

bool file_exists(const std::string &fname) { return !std::ifstream{fname}.fail(); }

hid_t hdf5_open_file(const std::string &fname, bool read_only) {
  hid_t id;
  if (read_only)
    id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  else
    id = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  return id;
}

hid_t hdf5_create_file(const std::string &fname) {
  return H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

hid_t hdf5_open_file(const std::string &fname, MPI_Comm communicator, bool read_only) {
  auto plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, communicator, MPI_INFO_NULL);
  hid_t id;
  if (read_only)
    id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, plist_id);
  else
    id = H5Fopen(fname.c_str(), H5F_ACC_RDWR, plist_id);
  H5Pclose(plist_id);
  return id;
}

hid_t hdf5_create_file(const std::string &fname, MPI_Comm communicator) {
  auto plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, communicator, MPI_INFO_NULL);
  auto id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
  return id;
}

bool hdf5_file_is_open(hid_t file_id) { return H5Fget_obj_count(file_id, H5F_OBJ_FILE) > 0; }

bool hdf5_dataset_exists(hid_t location, const std::string &dataset_name) {
  return H5Lexists(location, dataset_name.c_str(), H5P_DEFAULT) > 0;
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro
