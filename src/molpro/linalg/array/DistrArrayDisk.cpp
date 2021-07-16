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
    : DistrArray(distr->border().second, commun), m_distribution(std::move(distr)) {allocate_chunks();}

DistrArrayDisk::DistrArrayDisk() {allocate_chunks();}

DistrArrayDisk::DistrArrayDisk(const DistrArrayDisk& source)
    : DistrArray(source), m_distribution(std::make_unique<Distribution>(*source.m_distribution)) {
  DistrArray::copy(source);
  allocate_chunks();
}

DistrArrayDisk::DistrArrayDisk(DistrArrayDisk&& source) noexcept
    : DistrArray(source), m_distribution(std::move(source.m_distribution)) {
  using std::swap;
  if (source.m_allocated) {
    swap(m_allocated, source.m_allocated);
  }
  allocate_chunks();
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

void DistrArrayDisk::allocate_chunks(){
  this->chunks.emplace_back(chunk_size);
  this->chunks.emplace_back(chunk_size);
  std::future<void> dummy_future_1 = std::async(std::launch::async, []{  });
  std::future<void> dummy_future_2 = std::async(std::launch::async, []{  });
  this->chunk_futures.emplace_back(std::move(dummy_future_1));
  this->chunk_futures.emplace_back(std::move(dummy_future_2));
}

DistrArray::value_type DistrArrayDisk::dot(const DistrArray& y) const {
  std::cout << "Old dot\n";
  return dot(y);
}

void DistrArrayDisk::toggle_chunk(){
  this->curr_chunk = (this->curr_chunk + 1) % 2;
}

std::vector<DistrArray::value_type>& DistrArrayDisk::get_chunk(){
  int chunk_to_get = this->curr_chunk;
  this->chunk_futures[chunk_to_get].wait_for(std::chrono::seconds(0)) == std::future_status::ready;
  return this->chunks[chunk_to_get];
}

  // todo: this is loading both buffers when it should just be loading one
void DistrArrayDisk::load_buffers(int i){
  int buffer = this->curr_chunk;
  auto range = this->m_distribution->range(molpro::mpi::rank_global());
  size_t offset = range.first + i*this->chunk_size;
  int nextbuffer = (buffer + 1) % 2;
  const auto hi = std::min(offset + this->chunk_size, range.second);
  std::future<void> buffer_1, buffer_2;
  // load up the current buffer
  if (offset == range.first)
    get(offset, hi, this->chunks[buffer].data());
  else {
    //get(offset, hi, this->chunks[buffer].data()); // this but async
    this->chunk_futures[buffer] = std::async(std::launch::async, [this, offset, hi, buffer]{ this->get(offset, hi, this->chunks[buffer].data()); });
  }
  // load up the next buffer
  if (offset + this->chunk_size < range.second) {
    const auto hinext = std::min(offset + 2 * this->chunk_size, range.second);
    //get(offset + this->chunk_size, hinext, this->chunks[nextbuffer].data());
    this->chunk_futures[nextbuffer] = std::async(std::launch::async, [this, offset, hinext, nextbuffer]{ this->get( offset + this->chunk_size, hinext, this->chunks[nextbuffer].data() ); });
  }

}

bool DistrArrayDisk::in_chunks(int i){
  auto range = this->m_distribution->range(molpro::mpi::rank_global());
  return i*this->chunk_size + range.first < range.second;
}

// todo: this needs to be const, class variables need to be mutable (note: don't do this, make it own class instead)
// todo: buffer mgmt should be its own class that is instantiated here
// todo: do away with toggle_chunk(), instead toggle after load_buffers() (note: put toggle chunk in in_chunks)
// don't pass non-const references - these buffers don't belong to the user
// class can have an optional argument in constructor for chunk size

DistrArray::value_type DistrArrayDisk::dot(const DistrArray& y) {
  DistrArray::value_type result = 0;
  auto range = this->m_distribution->range(molpro::mpi::rank_global());
  for (size_t i = 0; in_chunks(i); i+=1) {
    size_t chunk_start = range.first + i*this->chunk_size;
    const auto chunk_end = std::min(chunk_start + this->chunk_size, range.second);
    DistrArrayDisk::load_buffers(i);
    result = std::inner_product(begin(get_chunk()),
                                begin(get_chunk()) + chunk_end - chunk_start,
                                y.local_buffer()->data() + chunk_start - range.first,
                                result);
    this->toggle_chunk();
  }
  return result;
}
DistrArray::value_type DistrArrayDisk::dot(const DistrArray::SparseArray& y) const { return DistrArray::dot(y); }

} // namespace molpro::linalg::array
