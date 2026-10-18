#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_H
#include <array>
#include <atomic>
#include <cstddef>
#include <future>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

#include <molpro/linalg/array/Span.h>
#include <molpro/mpi.h>

namespace molpro::linalg::array {

namespace util {

template <typename T>
struct is_std_array : std::false_type {};

template <typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

template <typename T>
constexpr bool is_std_array_v = is_std_array<T>::value;

template <typename T>
constexpr bool is_array_v = std::is_array<T>::value || is_std_array_v<T>;

template <typename T, typename = void>
struct is_allocatable : std::false_type {};

// We use the existence of a constructor taking a size-type as a proxy
// for whether the given type is "allocatable", that is it can be newly
// constructed by simply telling it how many elements it shall contain.
template <typename T>
struct is_allocatable<T, std::enable_if_t<std::is_constructible_v<T, std::size_t>>> : std::true_type {};

template <typename T>
constexpr bool is_allocatable_v = is_allocatable<T>::value;

} // namespace util

/*!
 * @brief Allocates an array of type T to hold the specified amount of items
 *
 * This function serves as a customization point for array types that don't support
 * allocation via a constructor that takes in a size. The default implementation will
 * throw an exception for such types and is hence only viable in case no such allocations
 * are required during any point of execution.
 *
 * @param size The number of objects the to-be-allocated array should encompass
 * @returns The allocated array
 *
 * @note There exists no hook for deallocation. It is assumed that the memory allocated by
 * specializations of this function is managed and tracked externally and freed once
 * this library has done all requested work (or in between solver invocations).
 */
template<typename T>
T allocate_array(std::size_t size) {
  if constexpr (util::is_allocatable_v<T>) {
    return T(size);
  }

  throw std::runtime_error("Required molpro::linalg::array::allocate_array for type '"
          + std::string(typeid(T).name())
          + "' - you have to provide an explicit template specialization for this function & type");
}

namespace util {
/*!
 * @brief Creates an array of type T able to hold size elements
 *
 * Fixed-size arrays are default constructed, types constructible from a size are constructed from it, and any other
 * type is created by the customization point allocate_array().
 */
template <typename T>
T allocate_like(std::size_t size) {
  if constexpr (is_array_v<T>)
    return T();
  else if constexpr (is_allocatable_v<T>)
    return T(size);
  else
    return allocate_array<T>(size);
}

//! A view of all the elements of a contiguous array
template <class R>
auto contiguous_span(R &r) {
  return Span<std::remove_pointer_t<decltype(std::data(r))>>(std::data(r), std::size(r));
}
} // namespace util

} // namespace molpro::linalg::array

namespace molpro::linalg::array::util {

//! Atomic lock allowing only one process to acquire it. Implemented using MPI3 RMA.
class LockMPI3 {
protected:
  MPI_Comm m_comm = MPI_COMM_NULL; //!< MPI communicator
  MPI_Win m_win = MPI_WIN_NULL;    //!< empty window handle
  bool m_locked = false;           //!< whether lock is active

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

//! Utility object that locks on creation and unlocks on destruction
class ScopeLock {
public:
  explicit ScopeLock(MPI_Comm comm) : lock{comm}, l{lock.scope()} {}

protected:
  LockMPI3 lock;
  decltype(std::declval<LockMPI3>().scope()) l;
};

template <typename Result = void>
class Task {
public:
  Task(std::future<Result> &&task) : m_task{std::move(task)} {};
  Task(Task &&other) = default;

  template <class Func, typename... Args>
  static Task create(Func &&f, Args &&...args) {
    return {std::async(std::launch::async, std::forward<Func>(f), std::forward<Args>(args)...)};
  }

  ~Task() {
    // Don't call wait() here: a moved-from Task has an invalid future, and
    // wait() throws future_error in that case. A throwing destructor in a
    // non-noexcept context would terminate the program.
    if (m_task.valid())
      m_task.wait();
  }

  //! Returns true if the task has completed
  bool test() { return m_task.wait_for(std::chrono::microseconds{1}) == std::future_status::ready; }
  //! wait for the task to complete and return its result
  Result wait() {
    if (!m_task.valid())
      throw std::future_error(std::future_errc::no_state);
    return m_task.get();
  }

protected:
  std::future<Result> m_task; //! Future holding the thread
};

} // namespace molpro::linalg::array::util

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_H
