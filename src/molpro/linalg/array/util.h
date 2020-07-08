#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
#include "molpro/linalg/array/DistrArray.h"
#include <memory>
#include <mpi.h>
#include <vector>

namespace molpro {
class Profiler;

namespace linalg {
namespace array {
namespace util {

//! Starts profiler on construction and stops it on destruction, if passed pointer is not null
struct ScopeProfiler {
  std::weak_ptr<molpro::Profiler> m_prof;
  std::string m_name;
  ScopeProfiler(const std::shared_ptr<molpro::Profiler> &prof, std::string name);
  ~ScopeProfiler();
  ScopeProfiler() = delete;
  ScopeProfiler(const ScopeProfiler &) = delete;
  ScopeProfiler &operator=(const ScopeProfiler &) = delete;
};

//! Atomic lock allowing only one process to acquire it. Implemented using MPI3 RMA.
class LockMPI3 {
protected:
  MPI_Comm m_comm = MPI_COMM_NULL; //! MPI communicator
  MPI_Win m_win = MPI_WIN_NULL;    //! empty window handle
  bool m_locked = false;           //!

public:
  //! Create the lock without acquiring it. This is collective and must be called by all processes in the communicator.
  LockMPI3(MPI_Comm comm);
  //! Release the lock and destroy it. This is collective and must be called by all processes in the communicator.
  ~LockMPI3();

  LockMPI3() = delete;
  LockMPI3(const LockMPI3 &) = delete;
  LockMPI3 &operator=(const LockMPI3 &) = delete;

  //! Acquire exclusive lock
  void lock();

  //! Release the lock
  void unlock();

protected:
  //! Proxy that locks on creation and unlocks on destruction. Useful for locking a scope.
  struct Proxy {
    LockMPI3 &m_lock;
    explicit Proxy(LockMPI3 &);
    Proxy() = delete;
    ~Proxy();

    bool m_deleted = false; //!< whether the lock was already deleted
  };

  //! Keep track of proxy object so that if lock is deleted, the proxy does not try to unlock.
  std::weak_ptr<Proxy> m_proxy;

public:
  //! Return a proxy object that acquires the lock on creation and releases it on destruction. Useful for locking a
  //! scope.
  std::shared_ptr<Proxy> scope();
};
} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
