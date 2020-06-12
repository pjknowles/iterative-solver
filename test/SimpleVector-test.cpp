#include "molpro/linalg/SimpleArray.h"
#include "test.h"
#include <cmath>

TEST(SimpleVector, copy_constructor) {
  using pv = molpro::linalg::SimpleArray<double>;
  pv v0(10001);
  for (size_t i = 0; i < v0.size(); i++) v0[i] = 2 * i + 1;
  bool result = true;
  for (size_t i = 0; i < 4; i++) {
//        std::cout <<"option "<<i<< ", copy-construct"<<std::endl;
    auto v2 = pv(v0, i);
//        std::cout <<"option "<<i<< ", copy-construct, reads="<<v2.m_cache.reads << " writes="<<v2.m_cache.writes<<std::endl;
    result &= v2 == v0;
//        std::cout <<"option "<<i<< ", compare, reads="<<v2.m_cache.reads << " writes="<<v2.m_cache.writes<<std::endl;
//        std::cout <<mpi_rank<<"source: "<<v0<<std::endl;
//        std::cout <<mpi_rank<<"copy: "<<v2<<std::endl;
//        std::cout <<mpi_rank<<"result: "<<result<<std::endl;
  }
//       std::cout <<"After Copy Constructor "<<SimpleVector<double>(v0,1)<<std::endl;
//       auto test = SimpleVector<double>(v0,1);
//       std::cout <<"After Copy Constructor "<<test<<std::endl;
  ASSERT_TRUE(result);
}

TEST(SimpleVector, dot_product) {
  molpro::linalg::SimpleArray<double> v0(10001);
  for (size_t i = 0; i < v0.size(); i++) v0[i] = 2 * i + 1;
//       std::cout << "v0=k<<v0<<std::endl;
  bool result = true;
  for (size_t i = 0; i < 4; i++) {
//       for (size_t i=2; i<3; i++) {
//       std::cout << "before copy to v1 v0="<<v0<<std::endl;
    auto v1 = molpro::linalg::SimpleArray<double>(v0, i);
//       std::cout << "after copy to v1 v0="<<v0<<std::endl;
//       std::cout << "after copy to v1 v1="<<v0<<std::endl;
//        for (auto j=i-i%2; j<=i; j++) {
    for (auto j = 0; j < 4; j++) {
//         for (auto j = 1; j < 2; j++) {
      auto v2 = molpro::linalg::SimpleArray<double>(v0, j);
//         std::cout << "v0="<<v0<<std::endl;
//         std::cout << "v2="<<v2<<std::endl;
//      std::cout << v0.size() << std::endl;
//      double l = v0.size();
//      double ans = l * (2 * l - 1) * (2 * l + 1)/3;
//      std::cout << i << "," << j << "v0.v0: " << v0.dot(v0) << ": " << ans << ": " << v0.dot(v0) - ans << std::endl;
//      std::cout << i << "," << j << "v1.v0: " << v1.dot(v0) << ": " << ans << ": " << v1.dot(v0) - ans << std::endl;
//      std::cout << i << "," << j << "v2.v0: " << v2.dot(v0) << ": " << ans << ": " << v2.dot(v0) - ans << std::endl;
//      std::cout << i << "," << j << "v2.v1: " << v2.dot(v1) << ": " << ans << ": " << v2.dot(v1) - ans << std::endl;
//      std::cout << i << "," << j << ": " << v2.dot(v1) << "=="
//                << v0.size() * (2 * v0.size() - 1) * (2 * v0.size() + 1) / 3 << std::endl;;
//      std::cout << i << "," << j << ": " << v2.dot(v2) << "==" << v0.size() * (2 * v0.size() - 1) * (2 * v0.size() + 1) / 3 << std::endl;;
//      result &= std::fabs(v2.dot(v2) - double(v0.size() * (2 * v0.size() - 1) * (2 * v0.size() + 1) / 3)) < 1e-1;
//      result &= std::fabs(v2.dot(v1) - double(v0.size() * (2 * v0.size() - 1) * (2 * v0.size() + 1) / 3)) < 1e-1;
      result &= v2.dot(v2) == v0.size() * (2 * v0.size() - 1) * (2 * v0.size() + 1) / 3;
      result &= v2.dot(v1) == v0.size() * (2 * v0.size() - 1) * (2 * v0.size() + 1) / 3;
    }
  }
  ASSERT_TRUE(result);
}

#ifndef none
TEST(SimpleVector, axpy) {
  molpro::linalg::SimpleArray<double> v0(10001);
  for (size_t i = 0; i < v0.size(); i++) v0[i] = 2 * i + 1;
  bool result = true;
  for (auto i = 0; i < 4; i++) {
    for (auto j = 0; j < 4; j++) {
// for (auto i = 1; i < 2; i++) {
//  for (auto j = 2; j < 3 ; j++) {
      auto v1 = molpro::linalg::SimpleArray<double>(v0, i);
      auto v2 = molpro::linalg::SimpleArray<double>(v0, j);
//         std::cout <<mpi_rank<< "v0="<<v0<<std::endl;
      v2.axpy(2, v1);
//   std::cout <<mpi_rank<< "after first axpy v2="<<v2<<std::endl;
      v2.axpy(-3, v1);
//   std::cout <<mpi_rank<< "after second axpy v2="<<v2<<std::endl;
//   std::cout <<mpi_rank<<"axpy: "<<i<<","<<j<<": "<< v2.dot(v2) <<"==0"<<std::endl;;
//         std::cout << v2.dot(v2) <<std::endl;
      result &= v2.dot(v2) < 1e-20;
//         std::cout << "result "<<result<<v2.dot(v2)<<std::endl;
//       std::cout << "reads="<<v2.m_cache.reads << " writes="<<v2.m_cache.writes<<std::endl;
    }
  }
  ASSERT_TRUE(result);
}

#endif
#ifndef none
TEST(SimpleVector, scal) {
  molpro::linalg::SimpleArray<double> v0(10000);
  for (size_t i = 0; i < v0.size(); i++) v0[i] = 2 * i + 1;
  bool result = true;
  for (auto i = 0; i < 4; i++) {
    auto v1 = molpro::linalg::SimpleArray<double>(v0, i);
    auto v2 = molpro::linalg::SimpleArray<double>(v0, i);
    v2.scal(3);
    v2.axpy(-3, v1);
//         std::cout << v2.dot(v2) <<std::endl;
    result &= v2.dot(v2) < 1e-20;
  }
  ASSERT_TRUE(result);
}

TEST(SimpleVector, zero) {
  molpro::linalg::SimpleArray<double> v0(10001);
  for (size_t i = 0; i < v0.size(); i++) v0[i] = 2 * i + 1;
  bool result = true;
  for (auto i = 0; i < 4; i++) {
    auto v1 = molpro::linalg::SimpleArray<double>(v0, i);
    auto v2 = molpro::linalg::SimpleArray<double>(v0, i);
    v2.scal(0);
//         std::cout << v2.dot(&v2) <<std::endl;
    result &= v2.dot(v2) < 1e-20;
  }
  ASSERT_TRUE(result);
}

TEST(SimpleVector,put) {
  for (auto i = 0; i < 4; i++) {
//  MPI_Barrier(MPI_COMM_WORLD);
//  usleep(100);
//  std::cout << "option "<<i<<std::endl;

molpro::linalg::SimpleArray<double> v0(10001, i);
    size_t segment_length = (v0.size() - 1) / mpi_size + 1;
    v0.scal(0);
#ifndef none
    double val = 99;
    size_t offset = v0.size() / 2;
//       std::cout << "after zero, v0: "<<v0<<std::endl;

//    std::cout <<mpi_rank<< "before put, v0: "<<v0<<std::endl;
//  if (offset>=segment_length*mpi_rank && offset < segment_length*(mpi_rank+1))
//  if (mpi_rank==3)
    v0.put(&val, 1, offset);
//       std::cout <<mpi_rank<< "after put, v0: "<<v0<<std::endl;
#ifndef none
    double val2 = 77;
//       std::cout <<mpi_rank<< "before get "<<std::endl;
    v0.get(&val2, 1, offset);
    if (offset >= segment_length * mpi_rank && offset < segment_length * (mpi_rank + 1)) {
//       std::cout <<mpi_rank<< "after get, v0: "<<v0<<std::endl;
//       std::cout <<mpi_rank<<"val2: "<<val2<<"=="<<val<<std::endl;
      ASSERT_TRUE(val == val2);
    }
#endif
#ifndef none
    auto test = v0.dot(v0);
//       std::cout <<mpi_rank<<"test"<<test<<"=="<<val*val<<std::endl;
    if (offset >= segment_length * mpi_rank && offset < segment_length * (mpi_rank + 1)) {
      ASSERT_TRUE(test == val * val);
    }
#endif
//   MPI_Barrier(MPI_COMM_WORLD); std::cout << mpi_rank<<"first test done"<<std::endl;
#endif

#ifndef none
    std::vector<double> w;
    for (size_t k = 0; k < w.size(); k++) w.push_back(k);
    v0.put(w.data(), w.size(), 0);
//       std::cout << "v0 "<<v0<<std::endl;
    std::vector<double> w1;
    w1.assign(w.size(), 0);
// v0.get(w1.data(),w1.mpi_size()/2,w1.mpi_size()/2); // move the cache window
    v0.get(w1.data(), w1.size(), 0);
//   for (auto i=0; i<w1.mpi_size(); i++) std::cout << "w1: "<<w1[i]<<std::endl;
    auto test2 = 0;
    for (auto k = (i > 1 ? segment_length * mpi_rank : 0);
         k < (w1.size() && (i > 1 ? (k < segment_length * (mpi_rank + 1)) : true)); k++)
      test2 += std::fabs(w[k] - w1[k]);
    ASSERT_TRUE(test2 < 1e-15);
#endif
  }

}
#endif
