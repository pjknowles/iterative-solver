#include "temp_file.h"
#include <fstream>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

std::string temp_file_name(const std::string& base_name, const std::string& suffix) {
  auto fname = base_name + suffix;
  size_t i = 0;
  while (!std::ifstream{fname}.fail()) {
    fname = base_name + "_" + std::to_string(i++) + suffix;
  }
  return fname;
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro