#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_DISKARRAY_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_DISKARRAY_H

#include "molpro/gci/array/DistrArray.h"

#include <hdf5>

namespace molpro {
namespace gci {
namespace array {

//! Copy of DistrArray stored on disk
class DiskArray {
public:
  DiskArray() = delete;
  DiskArray(const DiskArray &) = delete;
  DiskArray(DiskArray &&) = delete;
  DiskArray &operator=(const DiskArray &) = delete;
  DiskArray &operator=(DiskArray &&) = delete;

  DiskArray(hid_t hid);
  DiskArray(const std::string& file_name);
  DiskArray(ArrayBase &, std::string file_name);
  virtual ~DiskArray();

protected:
  hid_t m_hid;
  bool m_file_owner; //! whether instance created the file
};

DiskArray &operator<<(const DiskArray &, const Array &);
DiskArray &operator>>(const DiskArray &, const Array &);

/*!
 * \code{.cpp}
 *   auto x = Array(n, comm);
 *   initialize(x);
 *   auto d = DiskArray(hid);
 *   d << x;
 *   auto y = Array(n, comm);
 *   d >> y;
 *   auto solver = IterativeSolver(options); // options should specify location of temporary files
 *   //within iterative solver
 *   onDisk.emplace_back(x, this->next_temp)
 * \endcode
 */

} // namespace array
} // namespace gci
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_DISKARRAY_H
