#include <gtest/gtest.h>

#include <filesystem>
#include <numeric>
#include <string>
#include <vector>

#include "parallel_util.h"
#include <molpro/linalg/array/ArrayHandlerDDiskIterable.h>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util/gemm.h>
#include <molpro/mpi.h>

using molpro::linalg::array::ArrayHandlerDDiskIterable;
using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::Span;
using molpro::linalg::test::mpi_comm;

namespace {
//! Number of files in directory
size_t number_of_files(const std::filesystem::path& directory) {
  return std::distance(std::filesystem::directory_iterator(directory), std::filesystem::directory_iterator{});
}

//! Whether two communicators hold the same processes in the same order
bool same_processes(MPI_Comm a, MPI_Comm b) {
  int comparison;
  MPI_Comm_compare(a, b, &comparison);
  return comparison == MPI_IDENT or comparison == MPI_CONGRUENT;
}
} // namespace

/*!
 * New disk arrays made by ArrayHandlerDDiskIterable or by DistrArrayFile::CreateTempCopy() have their files in the
 * directory given, which is shared by the processes, and leave no files behind
 */
class TempCopyDirectoryF : public ::testing::Test {
public:
  static constexpr size_t dim = 31;
  TempCopyDirectoryF() : directory(std::filesystem::absolute("TempCopyDirectory-test")), values(dim) {
    std::iota(values.begin(), values.end(), 1.);
    MPI_Comm_rank(mpi_comm, &rank);
    if (rank == 0) {
      std::filesystem::remove_all(directory);
      std::filesystem::create_directories(directory);
    }
    MPI_Barrier(mpi_comm);
  }
  ~TempCopyDirectoryF() override {
    MPI_Barrier(mpi_comm);
    if (rank == 0)
      std::filesystem::remove_all(directory);
    MPI_Barrier(mpi_comm);
  }

  /*!
   * Checks that a has the content values, and that its file is in directory. DistrArrayFile unlinks its file as soon as
   * it is open, so the file itself is not seen in the directory.
   */
  template <class A>
  void check(const A& a) {
    const auto [lo, hi] = molpro::linalg::array::util::local_range(a);
    EXPECT_EQ(a.vec(), std::vector<double>(values.begin() + lo, values.begin() + hi));
  }
  //! Checks that no files are left in directory
  void check_removed() {
    MPI_Barrier(mpi_comm);
    EXPECT_EQ(number_of_files(directory), 0);
  }

  const std::filesystem::path directory;
  std::vector<double> values;
  int rank = 0;
};

TEST_F(TempCopyDirectoryF, DistrArrayFile_CreateTempCopy_contiguous) {
  {
    auto a = DistrArrayFile::CreateTempCopy(Span<const double>(values.data(), values.size()), mpi_comm,
                                            directory.string());
    EXPECT_TRUE(same_processes(a.communicator(), mpi_comm));
    EXPECT_EQ(a.size(), dim);
    check(a);
  }
  check_removed();
}

TEST_F(TempCopyDirectoryF, ArrayHandlerDDiskIterable) {
  {
    auto handler = ArrayHandlerDDiskIterable<DistrArrayFile, std::vector<double>>(directory.string());
    auto a = handler.copy(values);
    EXPECT_TRUE(same_processes(a.communicator(), molpro::mpi::comm_self()));
    check(a);
  }
  check_removed();
}
