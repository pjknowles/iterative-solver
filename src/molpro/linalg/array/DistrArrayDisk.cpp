#include "DistrArrayDisk.h"
#include "util.h"
#include "util/Distribution.h"

namespace molpro::linalg::array {
using util::Task;
namespace {

int mpi_rank(MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  return rank;
}

} // namespace

DistrArrayDisk::DistrArrayDisk(std::unique_ptr<Distribution> distr, MPI_Comm commun)
    : DistrArray(distr->border().second, commun), m_distribution(std::move(distr)) {}

DistrArrayDisk::DistrArrayDisk() = default;

DistrArrayDisk::DistrArrayDisk(const DistrArrayDisk& source)
    : DistrArray(source),
      m_distribution(source.m_distribution ? std::make_unique<Distribution>(*source.m_distribution) : nullptr) {
  if (source.m_allocated) {
    DistrArray::copy(source);
  }
}

DistrArrayDisk::DistrArrayDisk(DistrArrayDisk&& source) noexcept
    : DistrArray(source), m_distribution(std::move(source.m_distribution)) {
  using std::swap;
  if (source.m_allocated) {
    swap(m_allocated, source.m_allocated);
    swap(m_view_buffer, source.m_view_buffer);
    swap(m_owned_buffer, source.m_owned_buffer);
  }
}

DistrArrayDisk::~DistrArrayDisk() = default;

DistrArrayDisk::LocalBufferDisk::LocalBufferDisk(DistrArrayDisk& source) : m_source{source} {
  int rank = mpi_rank(source.communicator());
  index_type hi;
  std::tie(m_start, hi) = source.distribution().range(rank);
  m_size = hi - m_start;
  if (!source.m_allocated) {
    m_snapshot_buffer.resize(m_size);
    m_buffer = &m_snapshot_buffer[0];
    m_dump = true;
    source.get(start(), start() + size(), m_buffer);
  } else {
    m_buffer = source.m_view_buffer.data();
    m_dump = false;
  }
}

bool DistrArrayDisk::LocalBufferDisk::is_snapshot() { return !m_snapshot_buffer.empty(); }

DistrArrayDisk::LocalBufferDisk::~LocalBufferDisk() {
  if (m_dump)
    m_source.put(start(), start() + size(), m_buffer);
}

void DistrArrayDisk::flush() {
  if (!m_allocated)
    return;
  auto rank = mpi_rank(communicator());
  index_type lo, hi;
  std::tie(lo, hi) = distribution().range(rank);
  put(lo, hi, m_view_buffer.data());
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

} // namespace molpro::linalg::array
