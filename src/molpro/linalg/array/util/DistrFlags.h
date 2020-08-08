#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRFLAGS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRFLAGS_H
#include <mpi.h>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

/*!
 * @brief Distributed array of integer flags with one flag per process and a locking mechanism to ensure atomic access.
 */
class DistrFlags {
protected:
  struct Proxy;
  MPI_Comm m_comm = {};                    //!< mpi communicator
  MPI_Win m_win = MPI_WIN_NULL;            //! empty window handle
  int *m_data = nullptr;                   //!< window buffer with one element
  std::shared_ptr<bool> m_alive = nullptr; //!< flags whether the object is alive (true) or was destroyed (false)
public:
  DistrFlags() = default;
  DistrFlags(const DistrFlags &);
  DistrFlags(DistrFlags &&);
  DistrFlags &operator=(const DistrFlags &);
  DistrFlags &operator=(DistrFlags &&);
  DistrFlags(MPI_Comm comm) : m_comm(comm) {}
  ~DistrFlags();

  friend void swap(DistrFlags &x, DistrFlags &y);

  //! Whether distributed array was allocated.
  bool empty() const { return m_win == MPI_WIN_NULL; }

  //! Access flag on current process through a proxy, locking access for all other processes until the proxy is
  //! destroyed
  const Proxy access() const;
  //! Access flag on current process rank a proxy, locking access for all other processes until the proxy is destroyed
  const Proxy access(int rank) const;
  Proxy access();
  Proxy access(int rank);

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
    Proxy(DistrFlags &flags, int rank, std::shared_ptr<bool> status)
        : m_dflags{flags}, m_rank{rank}, m_dflags_alive{std::move(status)} {}
    //! value of the flag
    int value() const;

    //! set value of flag
    void set(int val);

    //! Which rank the proxy has access to
    int rank() const { return m_rank; }

  protected:
    //! Checks that the DistrFlags object is valid and raises logic_error if it is not
    void validate() const;
    DistrFlags &m_dflags;                 //!< reference to DistrFlags object
    int m_rank;                           //!< rank of process to lock
    std::shared_ptr<bool> m_dflags_alive; //!< status of m_dflags object, whether it is alive or was destroyed
  };
};
} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRFLAGS_H
