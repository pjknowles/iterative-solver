#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYHDF5_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYHDF5_H
#include "molpro/gci/array/DistrArray.h"
#include <hdf5.h>

namespace molpro {
namespace gci {
namespace array {

/*!
 * @brief Distributed array storing the buffer on disk using hdf5.
 *
 * Before using operations that require access to the buffer, the hdf5 file must be opened.
 * To avoid corruption on accidental termination the file must be closed when the data is no longer needed.
 */
class DistrArrayHDF5 : public DistrArray {
protected:
  //! distribution of array buffer among processes. Stores start index and size for each
  std::unique_ptr<util::Distribution> m_distribution;
  bool m_allocated = false; //!< whether file has been opened and storage reserved
  bool m_open = false;      //!< whether the file has been opened
  hid_t m_hid;              //!< hdf5 file handle
public:
  DistrArrayHDF5() = delete;
  DistrArrayHDF5(DistrArrayHDF5 &&other) = delete;
  DistrArrayHDF5 &operator=(const DistrArrayHDF5 &&) = delete;

  //! Save the buffer to file fname, overwriting it if it exists and overwrite is true
  DistrArrayHDF5(const std::string &fname, size_t dimension, MPI_Comm commun, bool overwrite=false, std::shared_ptr<Profiler> prof = nullptr);
  //! Assume that
  DistrArrayHDF5(hid_t hid, size_t dimension, MPI_Comm commun, std::shared_ptr<Profiler> prof = nullptr);
  //! Copy constructor allocates the buffer if source is not empty
  DistrArrayHDF5(const DistrArray &source);
  DistrArrayHDF5(const DistrArrayHDF5 &source);
  DistrArrayHDF5 &operator=(const DistrArray &source);
  DistrArrayHDF5 &operator=(const DistrArrayHDF5 &source);
  ~DistrArrayHDF5() override;
};

} // namespace array
} // namespace gci
} // namespace molpro
#endif // GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYHDF5_H
