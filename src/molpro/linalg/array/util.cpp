#include "util.h"

namespace molpro::linalg::array::util {

LockMPI3::LockMPI3(MPI_Comm comm) : m_comm{comm} {
  char* base = nullptr;
  MPI_Win_allocate(0, 1, MPI_INFO_NULL, m_comm, &base, &m_win);
}

LockMPI3::~LockMPI3() {
  if (!m_proxy.expired())
    m_proxy.lock()->m_deleted = true;
  unlock();
  MPI_Win_free(&m_win);
}

void LockMPI3::lock() {
  if (!m_locked) {
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, m_win);
    m_locked = true;
  }
}

void LockMPI3::unlock() {
  if (m_locked) {
    MPI_Win_unlock(0, m_win);
    m_locked = false;
  }
}

std::shared_ptr<LockMPI3::Proxy> LockMPI3::scope() {
  auto p = std::shared_ptr<LockMPI3::Proxy>{};
  if (m_proxy.expired()) {
    p = std::make_shared<Proxy>(*this);
    m_proxy = p;
  }
  return m_proxy.lock();
}

LockMPI3::Proxy::Proxy(LockMPI3& source) : m_lock(source) { m_lock.lock(); }

LockMPI3::Proxy::~Proxy() {
  if (!m_deleted)
    m_lock.unlock();
}
} // namespace molpro::linalg::array::util
