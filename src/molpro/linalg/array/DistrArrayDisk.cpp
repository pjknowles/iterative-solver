#include "DistrArrayDisk.h"
#include "util.h"
#include "util/Distribution.h"

namespace molpro {
namespace linalg {
namespace array {
using util::Task;
namespace {

int mpi_rank(MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  return rank;
}

} // namespace

DistrArrayDisk::LocalBufferDisk::LocalBufferDisk(DistrArrayDisk& source) : m_source{source} {
  int rank = mpi_rank(source.communicator());
  index_type hi;
  std::tie(m_start, hi) = source.distribution().range(rank);
  m_size = hi - m_start;
  if (!source.m_allocated) {
    m_snapshot_buffer.resize(m_size);
    m_buffer = &m_snapshot_buffer[0];
    m_dump = true;
    source.get(start(), start() + size() - 1, m_buffer);
  } else {
    m_buffer = source.m_buffer.data();
    m_dump = false;
  }
}

bool DistrArrayDisk::LocalBufferDisk::is_snapshot() { return !m_snapshot_buffer.empty(); }

DistrArrayDisk::LocalBufferDisk::~LocalBufferDisk() {
  if (m_dump)
    m_source.put(start(), start() + size() - 1, m_buffer);
}

void DistrArrayDisk::allocate_buffer() {
  if (m_allocated)
    return;
  auto rank = mpi_rank(communicator());
  index_type lo, hi;
  std::tie(lo, hi) = distribution().range(rank);
  size_t sz = hi - lo;
  if (m_allocated_buffer.size() < sz)
    m_allocated_buffer.resize(sz);
  m_buffer = Span<value_type>(&m_allocated_buffer[0], m_allocated_buffer.size());
  m_allocated = true;
}

void DistrArrayDisk::allocate_buffer(Span<value_type> buffer) {
  auto rank = mpi_rank(communicator());
  index_type lo, hi;
  std::tie(lo, hi) = distribution().range(rank);
  size_t sz = hi - lo;
  if (buffer.size() < sz)
    error("provided buffer is too small");
  if (m_allocated) {
    std::copy(begin(m_buffer), end(m_buffer), begin(buffer));
    free_buffer();
  }
  swap(m_buffer, buffer);
  m_allocated = true;
  if (!m_allocated_buffer.empty()) {
    m_allocated_buffer.clear();
    m_allocated_buffer.shrink_to_fit();
  }
}

void DistrArrayDisk::free_buffer() {
  m_buffer = Span<value_type>{};
  m_allocated_buffer.clear();
  m_allocated_buffer.shrink_to_fit();
  m_allocated = false;
}

void DistrArrayDisk::flush() {
  if (!m_allocated)
    return;
  auto rank = mpi_rank(communicator());
  index_type lo, hi;
  std::tie(lo, hi) = distribution().range(rank);
  put(lo, hi - 1, m_buffer.data());
}

bool DistrArrayDisk::empty() const { return m_allocated; }

std::unique_ptr<DistrArray::LocalBuffer> DistrArrayDisk::local_buffer() {
  return std::make_unique<LocalBufferDisk>(*this);
}

std::unique_ptr<const DistrArray::LocalBuffer> DistrArrayDisk::local_buffer() const {
  auto l = std::make_unique<LocalBufferDisk>(*const_cast<DistrArrayDisk*>(this));
  l->dump() = false;
  return l;
}

std::unique_ptr<Task<void>> DistrArrayDisk::tput(index_type lo, index_type hi, const value_type* data) {
  return std::make_unique<Task<void>>(
      Task<void>::create(std::launch::async, [lo, hi, data, this]() { this->put(lo, hi, data); }));
}

std::unique_ptr<Task<void>> DistrArrayDisk::tget(index_type lo, index_type hi, value_type* buf) {
  return std::make_unique<Task<void>>(
      Task<void>::create(std::launch::async, [lo, hi, buf, this]() { this->get(lo, hi, buf); }));
}

std::unique_ptr<Task<std::vector<DistrArrayDisk::value_type>>> DistrArrayDisk::tget(index_type lo, index_type hi) {
  return std::make_unique<Task<std::vector<value_type>>>(
      Task<std::vector<value_type>>::create(std::launch::async, [lo, hi, this]() { return this->get(lo, hi); }));
}

std::unique_ptr<Task<DistrArrayDisk::value_type>> DistrArrayDisk::tat(index_type ind) const {
  return std::make_unique<Task<value_type>>(Task<value_type>::create(
      std::launch::async, [ ind, this ]() -> auto { return this->at(ind); }));
}

std::unique_ptr<Task<void>> DistrArrayDisk::tset(index_type ind, value_type val) {
  return std::make_unique<Task<void>>(Task<void>::create(
      std::launch::async, [ ind, val, this ]() -> auto { return this->set(ind, val); }));
}

std::unique_ptr<Task<void>> DistrArrayDisk::tacc(DistrArray::index_type lo, DistrArray::index_type hi,
                                                 const DistrArray::value_type* data) {
  return std::make_unique<Task<void>>(Task<void>::create(
      std::launch::async, [ lo, hi, data, this ]() -> auto { return this->acc(lo, hi, data); }));
}

std::unique_ptr<Task<std::vector<DistrArrayDisk::value_type>>>
DistrArrayDisk::tgather(const std::vector<index_type>& indices) const {
  return std::make_unique<Task<std::vector<value_type>>>(Task<std::vector<value_type>>::create(
      std::launch::async, [&indices, this ]() -> auto { return this->gather(indices); }));
}

std::unique_ptr<Task<void>> DistrArrayDisk::tscatter(const std::vector<index_type>& indices,
                                                     const std::vector<value_type>& data) {
  return std::make_unique<Task<void>>(Task<void>::create(
      std::launch::async, [&indices, &data, this ]() -> auto { return this->scatter(indices, data); }));
}

std::unique_ptr<Task<void>> DistrArrayDisk::tscatter_acc(std::vector<index_type>& indices,
                                                         const std::vector<value_type>& data) {
  return std::make_unique<Task<void>>(Task<void>::create(
      std::launch::async, [&indices, &data, this ]() -> auto { return this->scatter_acc(indices, data); }));
}

std::unique_ptr<Task<std::vector<DistrArrayDisk::value_type>>> DistrArrayDisk::tvec() const {
  return std::make_unique<Task<std::vector<value_type>>>(Task<std::vector<value_type>>::create(
      std::launch::async, [this]() -> auto { return this->vec(); }));
}

std::unique_ptr<Task<std::unique_ptr<DistrArray::LocalBuffer>>> DistrArrayDisk::tlocal_buffer() {
  return std::make_unique<Task<std::unique_ptr<DistrArray::LocalBuffer>>>(
      Task<std::unique_ptr<DistrArray::LocalBuffer>>::create(
          std::launch::async, [this]() -> auto { return this->local_buffer(); }));
}

std::unique_ptr<Task<std::unique_ptr<const DistrArray::LocalBuffer>>> DistrArrayDisk::tlocal_buffer() const {
  return std::make_unique<Task<std::unique_ptr<const DistrArray::LocalBuffer>>>(
      Task<std::unique_ptr<const DistrArray::LocalBuffer>>::create(
          std::launch::async, [this]() -> auto { return this->local_buffer(); }));
}

} // namespace array
} // namespace linalg
} // namespace molpro
