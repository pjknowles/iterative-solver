#ifndef ITERATIVESOLVER_TEST_H
#define ITERATIVESOLVER_TEST_H
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
static int mpi_rank = 0;
static int mpi_size = 1;
class MPIEnvironment : public ::testing::Environment {
 public:
  virtual void SetUp() {
#ifdef HAVE_MPI_H
    int mpiError = MPI_Init(nullptr, nullptr);
    ASSERT_FALSE(mpiError);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_rank == 0) std::cout << mpi_size << " MPI process" << (mpi_size > 1 ? "es " : "") << std::endl;
    ::testing::TestEventListeners& listeners =
    ::testing::UnitTest::GetInstance()->listeners();
    if (mpi_rank != 0)
      delete listeners.Release(listeners.default_result_printer());
#endif
  }
  virtual void TearDown() {
#ifdef HAVE_MPI_H
    int mpiError = MPI_Finalize();
    ASSERT_FALSE(mpiError);
#endif
  }
  virtual ~MPIEnvironment() {}
};

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
  return RUN_ALL_TESTS();
}
#endif //ITERATIVESOLVER_TEST_H
