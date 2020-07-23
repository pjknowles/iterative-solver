#ifndef ITERATIVESOLVER_TEST_H
#define ITERATIVESOLVER_TEST_H
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
#ifdef HAVE_PPIDD_H
#include <ppidd.h>
#endif
static int mpi_rank = 0;
static int mpi_size = 1;
class MPIEnvironment : public ::testing::Environment {
public:
  virtual void SetUp() {
#ifdef HAVE_MPI_H
    int result;
    MPI_Initialized(&result);
#ifdef HAVE_PPIDD_H
    if (!result)
      PPIDD_Initialize(nullptr, nullptr, PPIDD_IMPL_DEFAULT);
    MPI_Initialized(&result);
    int mpiError = !result;
#else
    int mpiError = result ? 0 : MPI_Init(nullptr, nullptr);
#endif
    ASSERT_FALSE(mpiError);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_rank == 0)
      std::cout << mpi_size << " MPI process" << (mpi_size > 1 ? "es " : "") << std::endl;
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (mpi_rank != 0)
      delete listeners.Release(listeners.default_result_printer());
#endif
  }
  virtual void TearDown() {
#ifdef HAVE_MPI_H
#ifdef HAVE_PPIDD_H
    PPIDD_Finalize();
#else
    int mpiError = MPI_Finalize();
    ASSERT_FALSE(mpiError);
#endif
#endif
  }
  virtual ~MPIEnvironment() {}
};

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
  return RUN_ALL_TESTS();
}

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  for (const auto& s : v)
    o << " " << s;
  return o;
}

#endif // ITERATIVESOLVER_TEST_H
