#include "PHDF5Handle.h"

namespace molpro {
namespace linalg {
namespace array {
namespace util {
hid_t PHDF5Handle::_open_plist() {
  auto plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, m_comm, MPI_INFO_NULL);
  return plist_id;
}

PHDF5Handle::PHDF5Handle(MPI_Comm comm) : HDF5Handle{}, m_comm{comm} {}
PHDF5Handle::PHDF5Handle(std::string file, MPI_Comm comm) : HDF5Handle{std::move(file)}, m_comm{comm} {}
PHDF5Handle::PHDF5Handle(std::string file, std::string group, MPI_Comm comm)
    : HDF5Handle{std::move(file), std::move(group)}, m_comm{comm} {}
PHDF5Handle::PHDF5Handle(hid_t hid, MPI_Comm comm, bool transfer_ownership)
    : HDF5Handle(hid, transfer_ownership), m_comm{comm} {}
PHDF5Handle::PHDF5Handle(const PHDF5Handle& source, MPI_Comm comm) : PHDF5Handle{comm} { *this = source; }
PHDF5Handle::PHDF5Handle(const PHDF5Handle& source) : PHDF5Handle{source, source.m_comm} {}
PHDF5Handle& PHDF5Handle::operator=(const PHDF5Handle& source) {
  static_cast<HDF5Handle&>(*this) = static_cast<const HDF5Handle&>(source);
  return *this;
}
PHDF5Handle::PHDF5Handle(PHDF5Handle&& source, MPI_Comm comm) noexcept : PHDF5Handle{comm} {
  *this = std::forward<PHDF5Handle>(source);
}
PHDF5Handle::PHDF5Handle(PHDF5Handle&& source) noexcept : PHDF5Handle{source, source.m_comm} {}
PHDF5Handle& PHDF5Handle::operator=(PHDF5Handle&& source) noexcept {
  if (file_is_open())
    close_file();
  m_comm = source.m_comm;
  m_file_hid = source.m_file_hid;
  m_group_hid = source.m_group_hid;
  m_file_name = source.m_file_name;
  m_group_name = source.m_group_name;
  m_file_owner = source.m_file_owner;
  m_group_owner = source.m_group_owner;
  source.m_file_owner = false;
  source.m_group_owner = false;
  auto dummy = PHDF5Handle{source.m_comm};
  source = dummy;
  return *this;
}
PHDF5Handle::~PHDF5Handle() {
  HDF5Handle::close_file();
  if (m_erase_on_destroy) {
    int rank;
    MPI_Comm_rank(communicator(), &rank);
    if (rank == 0)
      if (file_exists(file_name()))
        std::remove(file_name().c_str());
    MPI_Barrier(communicator());
  }
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro
