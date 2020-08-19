#include "DistrFlags.h"
#include <stdexcept>

namespace {
int mpi_rank(MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  return rank;
}

} // namespace
namespace molpro {
namespace linalg {
namespace array {
namespace util {

DistrFlags::DistrFlags(MPI_Comm comm, int value) : m_comm{comm} {
  int size_of_type = sizeof(int);
  int *data = nullptr;
  MPI_Win_allocate(size_of_type, size_of_type, MPI_INFO_NULL, m_comm, &data, &m_win);
  *data = value;
  m_counter = std::make_shared<int>(0);
}

DistrFlags::~DistrFlags() {
  if (!empty()) {
    if (!m_counter || (m_counter && *m_counter != 0))
      MPI_Abort(m_comm, 1);
    MPI_Win_free(&m_win);
  }
}

DistrFlags::DistrFlags(const DistrFlags &source) : DistrFlags{source.m_comm} {
  auto value = source.access().get();
  access().replace(value);
}

DistrFlags::DistrFlags(DistrFlags &&source) noexcept : DistrFlags{} { swap(*this, source); }

DistrFlags &DistrFlags::operator=(const DistrFlags &source) {
  if (empty()) {
    DistrFlags t{source};
    swap(*this, t);
  } else {
    auto value = source.access().get();
    access().replace(value);
  }
  return *this;
}

DistrFlags &DistrFlags::operator=(DistrFlags &&source) noexcept {
  DistrFlags t{std::move(source)};
  swap(*this, t);
  return *this;
}

void swap(DistrFlags &x, DistrFlags &y) {
  using std::swap;
  swap(x.m_comm, y.m_comm);
  swap(x.m_win, y.m_win);
  swap(x.m_counter, y.m_counter);
}

bool DistrFlags::empty() const { return m_win == MPI_WIN_NULL; }

DistrFlags::Proxy DistrFlags::access(int rank) const {
  if (empty())
    MPI_Abort(m_comm, 2);
  return Proxy{m_comm, m_win, rank, m_counter};
}

DistrFlags::Proxy DistrFlags::access() const { return access(mpi_rank(m_comm)); }

DistrFlags::Proxy::Proxy(MPI_Comm comm, MPI_Win win, int rank, std::shared_ptr<int> counter)
    : m_comm{comm}, m_win{win}, m_rank{rank}, m_counter{std::move(counter)} {
  *m_counter += 1;
  if (*m_counter == 1)
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, m_rank, 0, m_win);
}

DistrFlags::Proxy::~Proxy() {
  if (*m_counter == 1)
    MPI_Win_unlock(m_rank, m_win);
  *m_counter -= 1;
}

int DistrFlags::Proxy::get() const {
  int val;
  MPI_Get(&val, 1, MPI_INT, m_rank, 0, 1, MPI_INT, m_win);
  return val;
}

int DistrFlags::Proxy::replace(int val) {
  int res;
  MPI_Fetch_and_op(&val, &res, MPI_INT, m_rank, 0, MPI_REPLACE, m_win);
  return res;
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro
