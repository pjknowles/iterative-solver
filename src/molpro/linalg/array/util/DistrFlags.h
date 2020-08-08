#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRFLAGS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRFLAGS_H
#include <memory>
#include <mpi.h>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

/*!
 * @brief Distributed array of integer flags with one flag per process and a locking mechanism to ensure atomic access.
 *
 * @note Copying creates a proxy to access flag values on their processor.
 *
 * @warning Access to the data of a given rank is done via a Proxy which locks the window buffer on that rank. Care must
 * be taken to ensure a combination of proxies does not lead to a deadlock, especially when mixed with copying.
 */
class DistrFlags {
protected:
  struct Proxy;
  MPI_Comm m_comm = {};         //!< mpi communicator
  MPI_Win m_win = MPI_WIN_NULL; //! empty window handle
  std::shared_ptr<int> m_counter =
      nullptr; //!< counter of Proxy objects created. They must all be destroyed before this can be destroyed
public:
  DistrFlags() = default;
  DistrFlags(const DistrFlags &source);
  DistrFlags(DistrFlags &&source) noexcept;
  DistrFlags &operator=(const DistrFlags &source);
  DistrFlags &operator=(DistrFlags &&source) noexcept;
  ~DistrFlags();

  //! Construct the distributed array with initial value of flag
  explicit DistrFlags(MPI_Comm comm, int value = 0);

  friend void swap(DistrFlags &x, DistrFlags &y);

  //! Whether distributed array was allocated.
  bool empty() const;

  //! Access flag on current process through a proxy, locking access for all other processes until the proxy is
  //! destroyed. Cannot be called on an empty object.
  Proxy access() const;
  //! Access flag on current process rank a proxy, locking access for all other processes until the proxy is destroyed.
  //! Cannot be called on an empty object.
  Proxy access(int rank) const;

protected:
  /*!
   * @brief Proxy for accessing the stored flags. The flag on a given rank is locked on construction and unlocked on
   * destruction.
   *
   * @warning The overlying DistrFlags object must be alive and allocated for Proxy to work. If it is destroyed or only
   * default constructed, than any access functions will throw a logic error.
   */
  class Proxy {
  public:
    /*!
     * @brief Construct Proxy using contents of an active DistrFlags object.
     *
     * Passing contents of DistrFlags object instead of a reference, allows movement of the overlying object without
     * invalidating the Proxy.
     *
     * @param comm communicator
     * @param win MPI window
     * @param rank rank of processor to access
     * @param counter number of proxies with this window
     */
    Proxy(MPI_Comm comm, MPI_Win win, int rank, std::shared_ptr<int> counter);
    ~Proxy();
    //! get value of the flag
    int get() const;

    //! replace value of flag @returns original value
    int replace(int val);

    //! Which rank the proxy has access to
    int rank() const { return m_rank; }

  protected:
    MPI_Comm m_comm = {};           //!< mpi communicator
    MPI_Win m_win = MPI_WIN_NULL;   //! empty window handle
    int m_rank;                     //!< rank of process to lock
    std::shared_ptr<int> m_counter; //!< counter of proxy objects created by overlying DistrFlags object.
  };
};
} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRFLAGS_H
