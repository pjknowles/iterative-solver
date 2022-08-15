#include "BufferManager.h"
#include "molpro/Profiler.h"
#include <molpro/linalg/array/util/Distribution.h>

namespace molpro::linalg::array::util {

BufferManager::BufferManager(const BufferManager::CVecRef<molpro::linalg::array::DistrArrayFile>& arrays,
                             size_t buffer_size, int number_of_buffers)
    : m_arrays(arrays), m_buffer_size(buffer_size), m_number_of_buffers(number_of_buffers),
      m_range(arrays.empty() ? std::pair<DistrArray::index_type, DistrArray::index_type>{0, 0}
                             : arrays.front().get().distribution().range(molpro::mpi::rank_global())),
      m_buffer(arrays.size() * number_of_buffers * buffer_size) {
  assert(number_of_buffers > 0 && number_of_buffers < 3);
}

Span<BufferManager::value_type> BufferManager::next(bool initial) {
  auto prof = molpro::Profiler::single()->push(std::string{"BufferManager::next("} + (initial ? "T" : "F") + ")");
  if (initial)
    m_current_segment = 0;
  else
    ++m_current_segment;
  const auto buffer_id = m_current_segment % m_number_of_buffers;
  const auto lo = m_range.first + m_current_segment * this->m_buffer_size;
  if (lo >= m_range.second)
    return Span<BufferManager::value_type>(nullptr, 0);
  const auto hi = std::min(lo + this->m_buffer_size, m_range.second);
  m_current_buffer_size = hi - lo;
  {
    if (m_number_of_buffers == 1 or lo == m_range.first) {
      for (size_t iarray = 0; iarray < m_arrays.size(); ++iarray) {
        assert(hi > lo);
        assert(buffer_id * m_buffer_size * m_arrays.size() + iarray * m_buffer_size + hi - lo < m_buffer.size());
        m_arrays[iarray].get().get(lo, hi,
                                   &m_buffer.at(buffer_id * m_buffer_size * m_arrays.size() + iarray * m_buffer_size));
      }
    } else {
      m_read_future.wait();
    }
  }

  const auto next_buffer_id = (m_current_segment + 1) % m_number_of_buffers;
  const auto next_lo = m_range.first + (m_current_segment + 1) * this->m_buffer_size;
  const auto next_hi = std::min(next_lo + this->m_buffer_size, m_range.second);
  if (m_number_of_buffers > 1 and next_lo < m_range.second) {
    m_read_future = std::async(std::launch::async, [this, next_lo, next_hi, next_buffer_id]() {
      for (size_t iarray = 0; iarray < m_arrays.size(); ++iarray)
        m_arrays[iarray].get().get(
            next_lo, next_hi, &m_buffer.at(next_buffer_id * m_buffer_size * m_arrays.size() + iarray * m_buffer_size));
    });
  }

  return Span<value_type>(&m_buffer.at(buffer_id * m_buffer_size * m_arrays.size()), m_buffer_size * m_arrays.size());
}

size_t BufferManager::buffer_stride() const { return m_buffer_size; }
size_t BufferManager::buffer_size() const { return m_current_buffer_size; }
size_t BufferManager::buffer_offset() const { return m_current_segment * m_buffer_size; }
} // namespace molpro::linalg::array::util
