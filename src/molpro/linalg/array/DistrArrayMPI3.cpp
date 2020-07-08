#include "DistrArrayMPI3.h"
#include "util.h"
#include <tuple>
namespace molpro::gci::array {

namespace {
int comm_size(MPI_Comm comm) {
  int res;
  MPI_Comm_size(comm, &res);
  return res;
}
} // namespace

DistrArrayMPI3::DistrArrayMPI3(size_t dimension, MPI_Comm commun, std::shared_ptr<Profiler> prof)
    : DistrArray(dimension, commun, std::move(prof)),
      m_distribution(std::make_unique<util::Distribution>(comm_size(m_communicator), m_dimension)) {}

DistrArrayMPI3::~DistrArrayMPI3() {
  if (m_allocated)
    free_buffer();
}

void DistrArrayMPI3::allocate_buffer() {
  if (!empty())
    return;
  int rank;
  MPI_Aint n;
  MPI_Comm_rank(m_communicator, &rank);
  std::tie(std::ignore, n) = m_distribution->proc_buffer[rank];
  MPI_Win_allocate(n, sizeof(value_type), MPI_INFO_NULL, m_communicator, m_base, m_win);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, *m_win);
}

bool DistrArrayMPI3::empty() const { return !m_allocated; }

DistrArrayMPI3::DistrArrayMPI3(const DistrArray& source)
    : DistrArray(source.size(), source.communicator(), nullptr),
      m_distribution(std::make_unique<util::Distribution>(comm_size(m_communicator), m_dimension)) {
  if (!source.empty()) {
    DistrArrayMPI3::allocate_buffer();
    DistrArray::copy(source);
  }
}

DistrArrayMPI3& DistrArrayMPI3::operator=(const DistrArray& source) {
  auto old_dim = m_dimension;
  m_dimension = source.size();
  m_communicator = source.communicator();
  m_distribution = std::make_unique<util::Distribution>(comm_size(m_communicator), m_dimension);
  if (!source.empty()) {
    if (empty())
      allocate_buffer();
    else if (old_dim != source.size()) {
      free_buffer();
      allocate_buffer();
    }
    copy(source);
  }
  return *this;
}
void DistrArrayMPI3::free_buffer() {
  if (!m_allocated)
    error("Attempting to free a buffer that was not allocated");
  MPI_Win_unlock_all(*m_win);
  MPI_Win_free(m_win);
  m_allocated = false;
}
} // namespace molpro::gci::array
