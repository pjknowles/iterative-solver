#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_PHDF5HANDLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_PHDF5HANDLE_H

#include "molpro/linalg/array/HDF5Handle.h"
#include <mpi.h>

namespace molpro {
namespace linalg {
namespace array {
namespace util {
/*!
 * @brief HDF5Handle that opens a file in parallel mode. Any operation that needs to open or close a file is collective.
 */
class PHDF5Handle : public HDF5Handle {
public:
  PHDF5Handle() = delete;
  //! @copydoc HDF5Handle::HDF5Handle()
  explicit PHDF5Handle(MPI_Comm comm);
  //! @copydoc HDF5Handle::HDF5Handle(std::string)
  PHDF5Handle(std::string file, MPI_Comm comm);
  //! @copydoc HDF5Handle::HDF5Handle(std::string, std::string)
  PHDF5Handle(std::string file, std::string group, MPI_Comm comm);
  //! @copydoc HDF5Handle::HDF5Handle(hid_t, bool)
  //! @param comm communicator
  PHDF5Handle(hid_t hid, MPI_Comm comm, bool transfer_ownership = false);
  ~PHDF5Handle() override;

  //! @copydoc HDF5Handle::HDF5Handle(const HDF5Handle&)
  //! @param comm assign a different communicator
  PHDF5Handle(const PHDF5Handle &source, MPI_Comm comm);
  //! @copydoc HDF5Handle::HDF5Handle(const HDF5Handle&)
  //! @note Uses the communicator from source
  PHDF5Handle(const PHDF5Handle &source);
  //! @copydoc HDF5Handle::operator=(const HDF5Handle&)
  //! @note The communicator is not copied.
  PHDF5Handle &operator=(const PHDF5Handle &source);
  //! @copydoc HDF5Handle::HDF5Handle(HDF5Handle&&)
  //! @param source source handle to move
  //! @param comm assigns a new communicator
  PHDF5Handle(PHDF5Handle &&source, MPI_Comm comm) noexcept;
  //! @copydoc HDF5Handle::HDF5Handle(HDF5Handle&&)
  //! @note Takes on communicator from source
  PHDF5Handle(PHDF5Handle &&source) noexcept;
  //! @copydoc HDF5Handle::operator=(HDF5Handle&&)
  //! @note Takes on communicator from source
  PHDF5Handle &operator=(PHDF5Handle &&source) noexcept;

  //! Communicator belonging to the handle. Communicator is bound on construction and cannot be changed.
  MPI_Comm communicator() const { return m_comm; }

protected:
  MPI_Comm m_comm = MPI_COMM_NULL; //!< processes in this communicator have parallel access to the file
  hid_t _open_plist() override;
};

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_PHDF5HANDLE_H
