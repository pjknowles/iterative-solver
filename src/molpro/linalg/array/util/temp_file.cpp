#include "temp_file.h"
#include "../util.h"
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <thread>
#include <unistd.h>

namespace molpro::linalg::array::util {
#ifdef HAVE_MPI_H
fs::path temp_file_name(const fs::path& base_name, const std::string& suffix, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  auto fname = fs::path{""}.native();
  {
#ifdef HAVE_MPI_H
    molpro::linalg::array::util::LockMPI3 lock(comm);
#endif
    fname = rank != 0 ? fs::path{""}.native() : temp_file_name(base_name, suffix).native();
  }
  int fnamesize = fname.size();
  MPI_Bcast(&fnamesize, 1, MPI_INT, 0, comm);
  fname.resize(fnamesize);
  MPI_Bcast((void*)fname.data(), fnamesize, MPI_CHAR, 0, comm);
  return fs::path{fname};
}
#endif

int s_temp_file_name_count = 0;
std::mutex s_mutex;
fs::path temp_file_name(const fs::path& base_name, const std::string& suffix) {
  auto lock = std::lock_guard(s_mutex);
  const std::string chars = "01234566789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const int length = 32;
  std::hash<std::thread::id> hasher;
  auto time = std::chrono::high_resolution_clock().now().time_since_epoch().count();
  auto tid = hasher(std::this_thread::get_id());
  auto pid = getpid();
  s_temp_file_name_count++;
  srand(int((time + tid + pid + s_temp_file_name_count) & std::numeric_limits<int>::max()));
  auto fname = base_name;
  for (int i = 0; i < length; i++)
    fname += chars[rand() % chars.size()];
  fname += suffix;
  return fname;
}

} // namespace molpro::linalg::array::util