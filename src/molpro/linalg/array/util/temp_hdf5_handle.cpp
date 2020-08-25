#include "temp_hdf5_handle.h"
#include "temp_file.h"
#include <cassert>

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

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro