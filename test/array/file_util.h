#ifndef LINEARALGEBRA_TEST_ARRAY_FILE_UTIL_H
#define LINEARALGEBRA_TEST_ARRAY_FILE_UTIL_H
#include <fstream>

namespace {
class GarbageCollector {
public:
  GarbageCollector() = default;
  GarbageCollector(const GarbageCollector&) = default;
  GarbageCollector(std::string fname) : file_name(std::move(fname)) {}
  ~GarbageCollector() { remove_test_files(); }
  void remove_test_files() {
    if (!file_name.empty())
      if (!std::ifstream{file_name}.fail())
        std::remove(file_name.c_str());
  }

  std::string file_name;
};
} // namespace

#endif // LINEARALGEBRA_TEST_ARRAY_FILE_UTIL_H
