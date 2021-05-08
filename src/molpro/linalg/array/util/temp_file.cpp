#include "temp_file.h"
#include <cstdlib>
#include <fstream>
#include <functional>
#include <thread>
#include <unistd.h>
#include <limits>

namespace molpro::linalg::array::util {
#ifdef HAVE_MPI_H
std::string temp_file_name(const std::string& base_name, const std::string& suffix, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  auto fname = rank != 0 ? "" : temp_file_name(base_name, suffix);
  int fnamesize = fname.size();
  MPI_Bcast(&fnamesize, 1, MPI_INT, 0, comm);
  fname.resize(fnamesize);
  MPI_Bcast((void*)fname.data(), fnamesize, MPI_CHAR, 0, comm);
  return fname;
}
#endif

std::string temp_file_name(const std::string& base_name, const std::string& suffix) {
  const std::string chars = "01234566789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const int length = 32;
  std::hash<std::thread::id> hasher;
  auto time = std::chrono::high_resolution_clock().now().time_since_epoch().count();
  auto tid = hasher(std::this_thread::get_id());
  auto pid = getpid();
  srand(int((time + tid + pid) & std::numeric_limits<int>::max()));
  auto fname = base_name;
  for (int i = 0; i < length; i++)
    fname += chars[rand() % chars.size()];
  fname += suffix;
  return fname;
}

} // namespace molpro::linalg::array::util