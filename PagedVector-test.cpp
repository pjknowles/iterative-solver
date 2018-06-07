#include "PagedVector.h"
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
   using namespace LinearAlgebra;
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI_H
  MPI_Init(&argc,&argv);
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if (rank==0) std::cout << size<<" MPI process"<<(size>1?"es ":"")<<std::endl;
//    std::cout << "MPI process number "<<rank<<std::endl;
  }
#else
 std::cout << "serial"<<std::endl;
#endif
  Catch::Session().run(argc,argv);
#ifdef HAVE_MPI_H
  MPI_Finalize();
#endif
  return 0;
}

      TEST_CASE("PagedVector copy constructor") {
    int rank;
#ifdef HAVE_MPI_H
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
     rank=0;
#endif
       PagedVector<double> v0(10001);
       for (size_t i=0; i<v0.size(); i++) v0[i]=2*i+1;
       bool result = true;
       for (size_t i=0; i<4; i++) {
//        std::cout <<"option "<<i<< ", copy-construct"<<std::endl;
        auto v2 = PagedVector<double>(v0,i);
//        std::cout <<"option "<<i<< ", copy-construct, reads="<<v2.m_cache.reads << " writes="<<v2.m_cache.writes<<std::endl;
        result &= v2==v0;
//        std::cout <<"option "<<i<< ", compare, reads="<<v2.m_cache.reads << " writes="<<v2.m_cache.writes<<std::endl;
//        std::cout <<rank<<"source: "<<v0<<std::endl;
//        std::cout <<rank<<"copy: "<<v2<<std::endl;
//        std::cout <<rank<<"result: "<<result<<std::endl;
       }
//       std::cout <<"After Copy Constructor "<<PagedVector<double>(v0,1)<<std::endl;
//       auto test = PagedVector<double>(v0,1);
//       std::cout <<"After Copy Constructor "<<test<<std::endl;
       REQUIRE(result);
      }

#ifndef none
      TEST_CASE("PagedVector dot product") {
       PagedVector<double> v0(10000);
       for (size_t i=0; i<v0.size(); i++) v0[i]=2*i+1;
//       std::cout << "v0="<<v0<<std::endl;
       bool result = true;
       for (size_t i=0; i<4; i++) {
//       std::cout << "before copy to v1 v0="<<v0<<std::endl;
         auto v1 = PagedVector<double>(v0,i);
//       std::cout << "after copy to v1 v0="<<v0<<std::endl;
//       std::cout << "after copy to v1 v1="<<v0<<std::endl;
        for (auto j=i-i%2; j<=i; j++) {
         auto v2 = PagedVector<double>(v0,j);
//         std::cout << "v0="<<v0<<std::endl;
//         std::cout << "v2="<<v2<<std::endl;
//         std::cout <<i<<","<<j<<": "<< v2.dot(v1) <<"=="<< v0.size()*(2*v0.size()-1)*(2*v0.size()+1)/3<<std::endl;;
         result &= v2.dot(v1) == v0.size()*(2*v0.size()-1)*(2*v0.size()+1)/3;
         result &= v2.dot(v2) == v0.size()*(2*v0.size()-1)*(2*v0.size()+1)/3;
        }
       }
       REQUIRE(result);
      }

      TEST_CASE("PagedVector axpy()") {
       PagedVector<double> v0(10000);
       for (size_t i=0; i<v0.size(); i++) v0[i]=2*i+1;
       bool result = true;
       for (auto i=0; i<4; i++) {
        auto v1 = PagedVector<double>(v0,i);
        auto v2 = PagedVector<double>(v0,i);
        v2.axpy(2,v1);
        v2.axpy(-3,v1);
//         std::cout << v2.dot(v2) <<std::endl;
         result &= v2.dot(v2) <1e-20;
//         std::cout << "result "<<result<<v2.dot(v2)<<std::endl;
//       std::cout << "reads="<<v2.m_cache.reads << " writes="<<v2.m_cache.writes<<std::endl;
        }
       REQUIRE(result);
      }

      TEST_CASE("PagedVector scal()") {
       PagedVector<double> v0(10000);
       for (size_t i=0; i<v0.size(); i++) v0[i]=2*i+1;
       bool result = true;
       for (auto i=0; i<4; i++) {
        auto v1 = PagedVector<double>(v0,i);
        auto v2 = PagedVector<double>(v0,i);
        v2.scal(3);
        v2.axpy(-3,v1);
//         std::cout << v2.dot(v2) <<std::endl;
         result &= v2.dot(v2) <1e-20;
        }
       REQUIRE(result);
      }

      TEST_CASE("PagedVector zero()") {
       PagedVector<double> v0(10000);
       for (size_t i=0; i<v0.size(); i++) v0[i]=2*i+1;
       bool result = true;
       for (auto i=0; i<4; i++) {
        auto v1 = PagedVector<double>(v0,i);
        auto v2 = PagedVector<double>(v0,i);
        v2.zero();
//         std::cout << v2.dot(&v2) <<std::endl;
         result &= v2.dot(v2) <1e-20;
        }
       REQUIRE(result);
      }

      TEST_CASE("PagedVector::put()") {
       PagedVector<double> v0(10000);
       v0.zero();
       double val=99;
       size_t offset=v0.size()/2;
//       std::cout << "after zero, v0: "<<v0<<std::endl;
       v0.put(&val,1,offset);
//       std::cout << "after put, v0: "<<v0<<std::endl;
       double val2=77;
//       std::cout << "before get "<<std::endl;
       v0.get(&val2,1,offset);
//       std::cout << "after get, v0: "<<v0<<std::endl;
//       std::cout <<"val2"<<val2<<std::endl;
       REQUIRE(val==val2);
       auto test = v0.dot(v0);
//       std::cout <<"test"<<test<<std::endl;
       REQUIRE(test==val*val);
      }
#endif

#ifdef none
      TEST_CASE("PagedVector::operator[]()") {
       PagedVector<double> v0(10);
       v0.zero();
       double val=99;
       size_t offset=1;
//       std::cout << "v0: "<<v0<<std::endl;
       v0.put(&val,1,offset);
//       std::cout << "v0: "<<v0<<std::endl;
//       std::cout <<v0[offset]<<std::endl;
//       auto spoiler=v0.dot(v0);
       auto test = v0[offset];
//       std::cout << "v0: "<<v0<<std::endl;
//       std::cout << test<<std::endl;
       REQUIRE(test==val);
      }
#endif
