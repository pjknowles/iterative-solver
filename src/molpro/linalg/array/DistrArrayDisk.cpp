#include "DistrArrayDisk.h"
#include "util.h"
#include "util/Distribution.h"
#include <iostream>

namespace molpro::linalg::array {
using util::Task;
namespace {

int mpi_rank(MPI_Comm comm) {
  if (comm == molpro::mpi::comm_global())
    return molpro::mpi::rank_global();
#ifdef HAVE_MPI_H
  int rank;
  MPI_Comm_rank(comm, &rank);
  return rank;
#endif
  throw std::logic_error("Attempt to access MPI communicator in serial mode");
}

} // namespace

DistrArrayDisk::DistrArrayDisk(std::unique_ptr<Distribution> distr, MPI_Comm commun)
    : DistrArray(distr->border().second, commun), m_distribution(std::move(distr)) {}

DistrArrayDisk::DistrArrayDisk() = default;

DistrArrayDisk::DistrArrayDisk(const DistrArrayDisk& source)
    : DistrArray(source), m_distribution(std::make_unique<Distribution>(*source.m_distribution)) {
  DistrArray::copy(source);
}

DistrArrayDisk::DistrArrayDisk(DistrArrayDisk&& source) noexcept
    : DistrArray(source), m_distribution(std::move(source.m_distribution)) {
  using std::swap;
  if (source.m_allocated) {
    swap(m_allocated, source.m_allocated);
  }
}

DistrArrayDisk::~DistrArrayDisk() = default;

DistrArrayDisk::LocalBufferDisk::LocalBufferDisk(DistrArrayDisk& source) : m_source{source} {
  std::cout << "LocalBufferDisk constructor " << this << std::endl;
  int rank = mpi_rank(source.communicator());
  index_type hi;
  std::tie(m_start, hi) = source.distribution().range(rank);
  m_size = hi - m_start;
  std::cout << "LocalBufferDisk " << this << " resizes from  " << m_snapshot_buffer.size() << " to " << m_size
            << std::endl;
  m_snapshot_buffer.resize(m_size);
  m_buffer = &m_snapshot_buffer[0];
  source.get(start(), start() + size(), m_buffer);
}

DistrArrayDisk::LocalBufferDisk::LocalBufferDisk(DistrArrayDisk& source, const Span<value_type>& buffer)
    : m_source{source} {
  int rank = mpi_rank(source.communicator());
  index_type hi;
  std::tie(m_start, hi) = source.distribution().range(rank);
  m_size = hi - m_start;
  if (m_size > buffer.size())
    source.error("LocalBufferDisk(): attempting to construct from a buffer that is too small");
  m_buffer = buffer.empty() ? nullptr : &buffer[0];
  source.get(start(), start() + size(), m_buffer);
}

DistrArrayDisk::LocalBufferDisk::~LocalBufferDisk() {
  std::cout << "LocalBufferDisk destructor " << this << ", size = " << m_snapshot_buffer.size() << std::endl;
  if (do_dump)
    m_source.put(start(), start() + size(), m_buffer);
}

const DistrArray::Distribution& DistrArrayDisk::distribution() const {
  if (!m_distribution)
    error("allocate buffer before asking for distribution");
  return *m_distribution;
}

std::unique_ptr<DistrArray::LocalBuffer> DistrArrayDisk::local_buffer() {
  return std::make_unique<LocalBufferDisk>(*this);
}

std::unique_ptr<const DistrArray::LocalBuffer> DistrArrayDisk::local_buffer() const {
  auto l = std::make_unique<LocalBufferDisk>(*const_cast<DistrArrayDisk*>(this));
  l->do_dump = false;
  return l;
}

std::unique_ptr<DistrArray::LocalBuffer> DistrArrayDisk::local_buffer(const Span<value_type>& buffer) {
  return std::make_unique<LocalBufferDisk>(*this, buffer);
}
std::unique_ptr<const DistrArray::LocalBuffer> DistrArrayDisk::local_buffer(const Span<value_type>& buffer) const {
  auto l = std::make_unique<LocalBufferDisk>(*const_cast<DistrArrayDisk*>(this), buffer);
  l->do_dump = false;
  return l;
}
DistrArray::value_type DistrArrayDisk::dot(const DistrArray& y) const {
  auto ylb = y.local_buffer();
  const size_t chunk_size = 3;
  DistrArray::value_type result{0};
  std::vector<std::vector<DistrArray::value_type>> thischunk;
  thischunk.emplace_back(chunk_size);
  thischunk.emplace_back(chunk_size);
  //  std::cout << "DistrArrayDisk::dot" << std::endl;
  auto range = this->m_distribution->range(molpro::mpi::rank_global());
  int buffer = 0;
  for (size_t offset = range.first; offset < range.second; offset += chunk_size) {
    int nextbuffer = (buffer + 1) % 2;
    const auto hi = std::min(offset + chunk_size, range.second);
    if (offset == range.first)
      get(offset, hi, thischunk[buffer].data());
    else {
      // complete asynchronous I/O
    }
    if (offset + chunk_size < range.second) {
      const auto hinext = std::min(offset + 2 * chunk_size, range.second);
      get(offset + chunk_size, hinext, thischunk[nextbuffer].data());
    }
    result = std::inner_product(begin(thischunk[buffer]), begin(thischunk[buffer]) + hi - offset,
                                ylb->data() + offset - range.first, result);
    buffer = nextbuffer;
  }
  return result;
}
DistrArray::value_type DistrArrayDisk::dot(const DistrArray::SparseArray& y) const { return DistrArray::dot(y); }

} // namespace molpro::linalg::array
