#include "DistrArrayDisk.h"
#include "util.h"
#include "util/Distribution.h"
#include <future>
#include <iostream>
#include <molpro/Profiler.h>

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

DistrArrayDisk::DistrArrayDisk() {}

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
//  std::cout << "LocalBufferDisk constructor " << this << std::endl;
  int rank = mpi_rank(source.communicator());
  index_type hi;
  std::tie(m_start, hi) = source.distribution().range(rank);
  m_size = hi - m_start;
//  std::cout << "LocalBufferDisk " << this << " resizes from  " << m_snapshot_buffer.size() << " to " << m_size
//            << std::endl;
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
//  std::cout << "LocalBufferDisk destructor " << this << ", size = " << m_snapshot_buffer.size() << std::endl;
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

DistrArray::value_type DistrArrayDisk::dot(const DistrArrayDisk& y) const{
  if (&y == this){
    throw std::invalid_argument("Cannot dot a DistrArrayDisk with itself");
  }
  return DistrArray::dot(y); //TODO: implement buffering in both DistrArrays
}

DistrArray::value_type DistrArrayDisk::dot(const DistrArray& y) const {
  if (&y == this){
    throw std::invalid_argument("Cannot dot a DistrArrayDisk with itself");
  }
  BufferManager buffer_manager = BufferManager(*this, this->m_buffer_size, BufferManager::Double);
  value_type result = 0;
  auto yy = y.local_buffer()->data();
  for (auto buffer = buffer_manager.begin(); buffer != buffer_manager.end(); yy += buffer->size(), ++buffer)
    result = std::inner_product(begin(*buffer), end(*buffer), yy, result);

#ifdef HAVE_MPI_H
  molpro::Profiler::single()->start("MPI_Allreduce");
  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, communicator());
  molpro::Profiler::single()->stop();
#endif
  return result;
}

DistrArray::value_type DistrArrayDisk::dot(const DistrArray::SparseArray& y) const { return DistrArray::dot(y); }

// BufferManager class functions

BufferManager::BufferManager(const DistrArrayDisk& distr_array_disk, size_t chunk_size,
                             BufferManager::buffertype buffers)
    : chunk_size(std::move(chunk_size)), distr_array_disk(distr_array_disk),
      range(distr_array_disk.distribution().range(molpro::mpi::rank_global())) {
  for (size_t buffer_count = 0; buffer_count < buffers; ++buffer_count)
    this->chunks.emplace_back(this->chunk_size);
}

Span<BufferManager::value_type> BufferManager::next(bool initial) {
  if (initial)
    curr_chunk = 0;

  const size_t offset = range.first + curr_chunk * this->chunk_size;
  const auto buffer_id = curr_chunk % chunks.size();
  std::vector<value_type>& buffer = chunks[buffer_id];
  if (offset >= range.second) return Span<BufferManager::value_type>(nullptr,0);
  if (chunks.size() == 1 or offset == range.first)
    this->distr_array_disk.get(offset, std::min(offset + this->chunk_size, range.second), buffer.data());
  else
    this->next_chunk_future.wait();

  const size_t next_offset = range.first + (curr_chunk + 1) * this->chunk_size;
  if (chunks.size() > 1 and next_offset < range.second) {
    const auto hi = std::min(next_offset + this->chunk_size, range.second);
    auto data = chunks[((curr_chunk + 1) % chunks.size())].data();
    this->next_chunk_future = std::async(
        std::launch::async, [this, next_offset, hi, data] { this->distr_array_disk.get(next_offset, hi, data); });
  }

  ++curr_chunk;
  return Span<value_type>(buffer.data(),
                          offset >= range.second ? 0 : std::min(size_t(chunk_size), range.second - offset));
}

} // namespace molpro::linalg::array
