#include "DistrArrayDisk.h"
namespace molpro {
namespace linalg {
namespace array {
namespace {

int mpi_rank(MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  return rank;
}

} // namespace

DistrArrayDisk::DistrArrayDisk(size_t dimension, MPI_Comm commun, std::shared_ptr<Profiler> prof)
    : DistrArray(dimension, commun, std::move(prof)) {}

DistrArrayDisk::LocalBufferDisk::LocalBufferDisk(DistrArrayDisk& source) : m_source{source} {
  if (!source.m_allocated) {
    m_start = m_size = 0;
    m_buffer = nullptr;
    return;
  }
  int rank;
  MPI_Comm_rank(source.communicator(), &rank);
  std::tie(m_start, m_size) = source.distribution().range(rank);
  bool get_data = false;
  if (!source.m_allocated) {
    m_snapshot_buffer.resize(m_size);
    m_buffer = m_snapshot_buffer.data();
    m_dump = true;
    get_data = true;
  } else {
    m_buffer = source.m_buffer.data();
    m_dump = false;
    get_data = !source.is_valid();
  }
  if (get_data) {
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, source.m_win);
    source.get(start(), start() + size() - 1, m_buffer);
    MPI_Win_unlock(rank, source.m_win);
  }
}

bool DistrArrayDisk::LocalBufferDisk::memory_view() { return m_snapshot_buffer.empty(); }

DistrArrayDisk::LocalBufferDisk::~LocalBufferDisk() {
  if (m_dump)
    m_source.put(start(), start() + size() - 1, m_buffer);
}

void DistrArrayDisk::allocate_buffer() {
  if (m_allocated)
    return;
  if (m_win == MPI_WIN_NULL) {
    int* base = nullptr;
    int size_of_type = sizeof(int);
    MPI_Win_allocate(size_of_type, size_of_type, MPI_INFO_NULL, m_communicator, &base, &m_win);
  }
  int rank;
  MPI_Comm_rank(communicator(), &rank);
  index_type _lo;
  size_t sz;
  std::tie(_lo, sz) = distribution().range(rank);
  if (m_allocated_buffer.size() < sz) {
    m_allocated_buffer.resize(sz);
    mark_validity(false);
  }
  m_buffer = Span<value_type>(m_allocated_buffer.data(), m_allocated_buffer.size());
  m_allocated = true;
}

void DistrArrayDisk::allocate_buffer(Span<value_type> buffer) {
  int rank;
  MPI_Comm_rank(communicator(), &rank);
  index_type _lo;
  size_t sz;
  std::tie(_lo, sz) = distribution().range(rank);
  if (buffer.size() < sz)
    error("provided buffer is too small");
  deallocate_buffer();
  m_buffer = buffer;
  m_allocated = true;
  mark_validity(false);
}

void DistrArrayDisk::deallocate_buffer() {
  m_buffer = Span<value_type>{};
  m_allocated_buffer.clear();
  m_allocated_buffer.shrink_to_fit();
  m_allocated = false;
  mark_validity(false);
}

void DistrArrayDisk::mark_validity(bool valid) {
  if (m_win == MPI_WIN_NULL)
    return;
  int buffer = valid ? 1 : 0;
  auto rank = mpi_rank(communicator());
  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, m_win);
  MPI_Put(&buffer, 1, MPI_INT, rank, 0, 1, MPI_INT, m_win);
  MPI_Win_unlock(rank, m_win);
}
bool DistrArrayDisk::is_valid() {
  if (m_win == MPI_WIN_NULL)
    return false;
  int buffer = 0;
  auto rank = mpi_rank(communicator());
  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, m_win);
  MPI_Get(&buffer, 1, MPI_INT, rank, 0, 1, MPI_INT, m_win);
  MPI_Win_unlock(rank, m_win);
  return buffer;
}

void DistrArrayDisk::flush() {
  if (!m_allocated)
    return;
  auto rank = mpi_rank(communicator());
  index_type _lo;
  size_t sz;
  std::tie(_lo, sz) = distribution().range(rank);
  if (is_valid())
    put(_lo, _lo + sz - 1, m_buffer.data());
  else
    get(_lo, _lo + sz - 1, m_buffer.data());
}

bool DistrArrayDisk::empty() const { return m_allocated; }

} // namespace array
} // namespace linalg
} // namespace molpro
