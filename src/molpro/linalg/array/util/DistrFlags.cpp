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
    if (m_counter || !m_counter && *m_counter != 0)
      MPI_Abort(m_comm, 1);
    MPI_Win_free(&m_win);
  }
}

DistrFlags::DistrFlags(const DistrFlags &source) : DistrFlags{source.m_comm} {
  auto value = source.access().value();
  access().set(value);
}

DistrFlags::DistrFlags(DistrFlags &&source) noexcept : DistrFlags{} { swap(*this, source); }

DistrFlags &DistrFlags::operator=(const DistrFlags &source) {
  if (empty()) {
    DistrFlags t{source};
    swap(*this, t);
  } else {
    auto value = source.access().value();
    access().set(value);
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

int DistrFlags::Proxy::value() const {}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro
