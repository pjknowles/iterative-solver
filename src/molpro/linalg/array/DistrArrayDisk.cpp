#include "DistrArrayDisk.h"
#include "util.h"
#include "util/Distribution.h"
#include <future>
#include <iostream>
#include <numeric>

namespace molpro::linalg::array {
using util::Task;
DistrArrayDisk::DistrArrayDisk(std::unique_ptr<Distribution> distr, MPI_Comm commun)
    : DistrArray(distr->border().second, commun), m_distribution(std::move(distr)) {}

DistrArrayDisk::DistrArrayDisk() = default;

DistrArrayDisk::DistrArrayDisk(const DistrArray& source)
    : DistrArray(source), m_distribution(std::make_unique<Distribution>(source.distribution())) {
  copy(source);
}

DistrArrayDisk::DistrArrayDisk(DistrArrayDisk&& source) noexcept
    : DistrArray(std::move(source)), m_distribution(std::move(source.m_distribution)) {
  using std::swap;
  if (source.m_allocated) {
    swap(m_allocated, source.m_allocated);
  }
}

DistrArrayDisk::~DistrArrayDisk() = default;

const DistrArray::Distribution& DistrArrayDisk::distribution() const {
  if (!m_distribution)
    error("allocate buffer before asking for distribution");
  return *m_distribution;
}

std::unique_ptr<DistrArray::LocalBuffer> DistrArrayDisk::local_buffer() {
  throw std::logic_error("DistrArrayDisk::local_buffer(): an array on disk has no local buffer");
}

std::unique_ptr<const DistrArray::LocalBuffer> DistrArrayDisk::local_buffer() const {
  throw std::logic_error("DistrArrayDisk::local_buffer(): an array on disk has no local buffer");
}

std::string DistrArrayDisk::temp_copy_directory(const DistrArray& source, const std::string& directory) {
  if (not directory.empty())
    return directory;
  if (const auto* disk = dynamic_cast<const DistrArrayDisk*>(&source))
    if (auto source_directory = disk->file_directory(); not source_directory.empty())
      return source_directory;
  return ".";
}

} // namespace molpro::linalg::array
