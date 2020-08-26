#include "temp_phdf5_handle.h"
#include "temp_file.h"
#include "temp_hdf5_handle.h"

#include <cassert>
#include <iostream>
namespace molpro {
namespace linalg {
namespace array {
namespace util {

PHDF5Handle temp_phdf5_handle(const std::string &base_name, MPI_Comm comm) {
  auto fname = temp_file_name(base_name, ".hdf5");
  auto t = PHDF5Handle(fname, comm);
  assert(t.set_erase_file_on_destroy(true));
  return t;
}

PHDF5Handle temp_phdf5_handle_group(const PHDF5Handle &handle, const std::string &base_name, MPI_Comm comm) {
  if (comm == MPI_COMM_NULL)
    comm = handle.communicator();
  std::string group_name;
  try {
    group_name = temp_group_name(handle, base_name);
  } catch (std::runtime_error &err) {
    std::cerr << "temp_phdf5_handle_group: " << err.what() << std::endl;
    MPI_Abort(comm, 1);
  }
  auto t = PHDF5Handle(handle);
  t.open_file(HDF5Handle::Access::read_write);
  t.open_group(group_name);
  if (not t.set_erase_group_on_destroy(true)) {
    std::cerr << "temp_phdf5_handle_group: failed to set_erase_group_on_destroy on the copy of handle" << std::endl;
    MPI_Abort(comm, 1);
  }
  return t;
}
} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro
