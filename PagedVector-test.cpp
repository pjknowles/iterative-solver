#include "IterativeSolver/PagedVector.h"
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
   using namespace LinearAlgebra;
#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if (rank==0) std::cout << size<<" MPI process"<<(size>1?"es ":"")<<std::endl;
//    std::cout << "MPI process number "<<rank<<std::endl;
  }
#endif
  Catch::Session().run(argc,argv);
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}

      TEST_CASE("PagedVector copy constructor") {
       PagedVector<double> v0(10000);
       for (auto i=0; i<v0.size(); i++) v0[i]=2*i+1;
       bool result = true;
       for (auto i=0; i<4; i++)
        result &= PagedVector<double>(v0,i)==v0;
       REQUIRE(result);
      }

      TEST_CASE("PagedVector dot product") {
       PagedVector<double> v0(10000);
       for (auto i=0; i<v0.size(); i++) v0[i]=2*i+1;
       bool result = true;
       for (auto i=0; i<4; i++) {
        auto v1 = PagedVector<double>(v0,i);
        for (auto j=i-i%2; j<=i; j++) {
         auto v2 = PagedVector<double>(v0,j);
//         std::cout <<i<<","<<j<<": "<< v2.dot(&v1) <<"=="<< v0.size()*(2*v0.size()-1)*(2*v0.size()+1)/3<<std::endl;;
         result &= v2.dot(&v1) == v0.size()*(2*v0.size()-1)*(2*v0.size()+1)/3;
        }
       }
       REQUIRE(result);
      }

//      TEST_CASE("PagedVector 0") {
//       REQUIRE(PagedVectorTest<double>(2,0).status);
//      }
//      TEST_CASE("PagedVector offline") {
//       REQUIRE(PagedVectorTest<double>(2,LINEARALGEBRA_CLONE_ADVISE_OFFLINE).status);
//      }
//      TEST_CASE("PagedVector distributed") {
//       REQUIRE(PagedVectorTest<double>(2,LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED).status);
//      }
//      TEST_CASE("PagedVector offline and distributed") {
//       REQUIRE(PagedVectorTest<double>(2,LINEARALGEBRA_CLONE_ADVISE_OFFLINE&&LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED).status);
//      }
