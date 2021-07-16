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
  auto range = this->m_distribution->range(molpro::mpi::rank_global());
  BufferManager buffer_manager = BufferManager(this, range);
  DistrArray::value_type result = 0;
  for (size_t i = 0; buffer_manager.in_chunks(i); i+=1) {
    size_t chunk_start = range.first + i*buffer_manager.chunk_size;
    //const auto chunk_end = std::min(chunk_start + buffer_manager.chunk_size, range.second);
    buffer_manager.load_buffers(i);
    result = std::inner_product(begin(buffer_manager.get_chunk()),
                                begin(buffer_manager.get_chunk()) + buffer_manager.get_buffer_end(i),
                                y.local_buffer()->data() + chunk_start - range.first,
                                result);
  }
  return result;
}

DistrArray::value_type DistrArrayDisk::dot(const DistrArray::SparseArray& y) const { return DistrArray::dot(y); }

// BufferManager class functions

BufferManager::BufferManager(const DistrArrayDisk *distr_array_disk, std::pair<size_t, size_t> range,
                             const size_t chunk_size)
  : chunk_size(chunk_size), distr_array_disk(distr_array_disk), range(range){
    this->allocate_chunks();
  }  

void BufferManager::allocate_chunks(){
  this->chunks.emplace_back(this->chunk_size);
  this->chunks.emplace_back(this->chunk_size);
  this->next_chunk_future = std::async(std::launch::async, []{  });
}

std::vector<DistrArray::value_type>& BufferManager::get_chunk(){
  return this->chunks[this->curr_chunk];
}

void BufferManager::load_buffers(int i){
  auto range = this->range;
  size_t offset = range.first + i*this->chunk_size;
  int next_chunk = (this->curr_chunk + 1) % 2;
  // load up the current buffer (first iteration only)
  if (offset == range.first){
    const auto hi = std::min(offset + this->chunk_size, range.second);
    this->distr_array_disk->get(offset, hi, this->chunks[this->curr_chunk].data());
  }
  // load up the next buffer
  if (offset + this->chunk_size < range.second) {
    const auto hinext = std::min(offset + 2 * this->chunk_size, range.second);
    this->next_chunk_future = std::async(std::launch::async, [this, offset, hinext, next_chunk]{
      this->distr_array_disk->get( offset + this->chunk_size, hinext, this->chunks[next_chunk].data() ); 
    });
  }

}

bool BufferManager::in_chunks(int i){
  auto range = this->range;
  this->curr_chunk = (this->curr_chunk + 1) % 2;
  this->next_chunk_future.wait_for(std::chrono::seconds(0));// == std::future_status::ready;
  return i*this->chunk_size + range.first < range.second;
}

size_t BufferManager::get_buffer_end(int i){
    size_t chunk_start = this->range.first + i*this->chunk_size;
    const auto chunk_end = std::min(chunk_start + this->chunk_size, this->range.second);
    return chunk_end - chunk_start;

}

} // namespace molpro::linalg::array
