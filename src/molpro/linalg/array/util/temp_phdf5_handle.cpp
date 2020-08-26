#include "temp_phdf5_handle.h"
#include "temp_file.h"

#include <cassert>
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

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro
