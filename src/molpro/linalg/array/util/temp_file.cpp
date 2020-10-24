#include "temp_file.h"
#include <fstream>
#include <functional>
#include <stdlib.h>
#include <thread>
#include <unistd.h>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

std::string temp_file_name(const std::string& base_name, const std::string& suffix) {
  const std::string chars = "01234566789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const int length = 32;
  std::hash<std::thread::id> hasher;
  auto time = std::chrono::high_resolution_clock().now().time_since_epoch().count();
  auto tid = hasher(std::this_thread::get_id());
  auto pid = getpid();
  srand(time + tid + pid);
  auto fname = base_name;
  for (int i = 0; i < length; i++)
    fname += chars[rand() % chars.size()];
  fname += suffix;
  return fname;
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro