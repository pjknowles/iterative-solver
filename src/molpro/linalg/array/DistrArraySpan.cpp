#include "DistrArraySpan.h"
#include "util/Distribution.h"
#include <string>
#include <memory>

namespace molpro::linalg::array {

namespace {
int mpi_size(MPI_Comm comm) {
  int rank;
  MPI_Comm_size(comm, &rank);
  return rank;
}
}

DistrArraySpan::DistrArraySpan(size_t dimension, MPI_Comm commun)
    : DistrArraySpan(std::make_unique<Distribution>(
      util::make_distribution_spread_remainder<index_type>(dimension, mpi_size(commun))), commun) {}
      
DistrArraySpan::DistrArraySpan(std::unique_ptr<Distribution> distribution, MPI_Comm commun)
    : DistrArray(distribution->border().second, commun), m_distribution(std::move(distribution)) {
  if (m_distribution->border().first != 0)
    DistrArray::error("Distribution of array must start from 0");
}

// DistrArraySpan::~DistrArraySpan() {}

DistrArraySpan::DistrArraySpan(const DistrArraySpan& source)
    : DistrArray(source.size(), source.communicator()),
       m_distribution(source.m_distribution ? std::make_unique<Distribution>(*source.m_distribution) : nullptr) {
  if (!source.empty()) {
    DistrArraySpan::allocate_buffer(*source.m_buffer); // TODO: is object slicing OK here??
  }
}

DistrArraySpan::DistrArraySpan(const DistrArray& source)
    : DistrArray(source), m_distribution(std::make_unique<Distribution>(source.distribution())) {
  if (!source.empty()) {
    DistrArraySpan::allocate_buffer(Span<value_type>(&(*source.local_buffer())[0], source.size()));
  }
}

DistrArraySpan::DistrArraySpan(DistrArraySpan&& source) noexcept
    : DistrArray(source.m_dimension, source.m_communicator), m_buffer(std::move(source.m_buffer)),
      m_allocated(source.m_allocated), m_distribution(std::move(source.m_distribution)) {
  source.m_buffer = nullptr;
  source.m_allocated = false;
}

DistrArraySpan& DistrArraySpan::operator=(const DistrArraySpan& source) {
  if (this == &source)
    return *this;
  if (source.empty() || empty() || !compatible(source)) {
    free_buffer();
    DistrArraySpan t{source};
    swap(*this, t);
  } else {
    allocate_buffer(*source.m_buffer);  // TODO: OK???
  }
  return *this;
}

DistrArraySpan& DistrArraySpan::operator=(DistrArraySpan&& source) noexcept {
  DistrArraySpan t{std::move(source)};
  swap(*this, t);
  return *this;
}

void swap(DistrArraySpan& a1, DistrArraySpan& a2) noexcept {
  using std::swap;
  swap(a1.m_dimension, a2.m_dimension);
  swap(a1.m_communicator, a2.m_communicator);
  swap(a1.m_distribution, a2.m_distribution);
  swap(a1.m_allocated, a2.m_allocated);
  swap(a1.m_buffer, a2.m_buffer);
}

void DistrArraySpan::free_buffer() {
  if (m_allocated) {
    m_buffer = nullptr;
    m_allocated = false;
  }
}

void DistrArraySpan::allocate_buffer() {} //TODO: Doesn't make sense to have it really...

void DistrArraySpan::allocate_buffer(Span<value_type> buffer) {
  if (!empty())
    return;
  if (!m_distribution)
    error("Cannot allocate an array without distribution");
  m_buffer = std::make_unique<LocalBufferSpan>(buffer);
  int rank = MPI_Comm_rank(m_communicator, &rank);
  index_type lo, hi;
  std::tie(lo, hi) = m_distribution->range(rank);
  size_t n = hi - lo;
  if (m_buffer->size() < n)
    error("Specified external buffer is too small");
  m_allocated = true;
}

DistrArraySpan::LocalBufferSpan::LocalBufferSpan(Span<value_type>& source) {
  // TODO: is there a smarter way of doing this (i.e. using Span copy-constructor)???
  m_buffer = source.data();
  m_size = source.size();
}

bool DistrArraySpan::empty() const { return !m_allocated; }

const DistrArray::Distribution& DistrArraySpan::distribution() const {
  if (!m_distribution)
    error("allocate buffer before asking for distribution");
  return *m_distribution;
}


} // molpro::linalg::array