#include "molpro/linalg/array/OutOfCoreArray.h"
#include "test.h"
#include <cmath>

TEST(OutOfCoreArray, copy_constructor) {
  using ocv = molpro::linalg::OutOfCoreArray<double, 100>;
  std::vector<double> v0(10001);
  for (size_t i = 0; i < v0.size(); i++) v0[i] = 2 * i + 1;
  bool result = true;
  ocv v1(v0.data(),v0.size()); // Constructor from external buffer
  ocv v2(v1); // Copy-constructor (writes on disk)
  result &= v2 == v1;
  //std::cout <<mpi_rank<<" source: ";
  //for (double i : v0) std::cout<<" "<<i;
  //std::cout<<std::endl;
  //std::cout <<mpi_rank<<" copy-1: "<<v1<<std::endl;
  //std::cout <<mpi_rank<<" copy-2: "<<v2<<std::endl;
  ASSERT_TRUE(result);
}

TEST(OutOfCoreArray, pass_through) { //Do we really need it?
  std::vector<double> v(10001);
  for (size_t i = 0; i < v.size(); ++i) v[i] = 2 * i + 1;
  auto v1 = molpro::linalg::OutOfCoreArray<double,100>(v.data(), v.size());
  auto& w = v1;
  double r = 0;
  bool ad = true;
  for (size_t i = 0; i < v.size(); ++i) {
    ad = ad && ((w[i]) == (v[i]));
    r += std::fabs(v[i] - w[i]);
  }
  ASSERT_TRUE(ad);
  ASSERT_TRUE(std::fabs(r) < 1e-15);
}

TEST(OutOfCoreArray, equal) {
    using ocv = molpro::linalg::OutOfCoreArray<double, 100>;
    std::vector<double> v(10001);
    for (size_t i = 0; i < v.size(); ++i) v[i] = 2 * i + 1;
    ocv v0(v.data(),v.size());
    ocv v1(v0);
    std::vector<double> t1(10001), t2(10001);
    for (size_t i = 0; i < t1.size(); ++i) t1[i] = i;
    for (size_t i = 0; i < t2.size(); ++i) t2[i] = i;
    ocv v2(t1.data(),t1.size());
    ocv v3(t2.data(),t2.size());
    bool result = true;
    result &= !(v0 == v2);
    v2 = v0;
    result &= v2.dot(v0) == v0.dot(v0);
    v3 = v1;
    result &= v3.dot(v1) == v1.dot(v1);
    ASSERT_TRUE(result);
}

TEST(OutOfCoreArray, dot) { //Also tests "==" operator
    using ocv = molpro::linalg::OutOfCoreArray<double, 100>;
    std::vector<double> v(10001);
    for (size_t i = 0; i < v.size(); ++i) v[i] = 2 * i + 1;
    ocv v0(v.data(),v.size());
    ocv v1(v0);
    bool result = true;
    result &= v0 == v1;
    /* Case of equal vectors */
    result &= v0.dot(v1) == v.size() * (2 * v.size() - 1) * (2 * v.size() + 1) / 3; //First in memory
    result &= v1.dot(v0) == v.size() * (2 * v.size() - 1) * (2 * v.size() + 1) / 3; //First on disk
    /* Case of different vectors */
    std::vector<double> t1(10001),t2(10001);
    for (size_t i = 0; i < t1.size(); ++i) t1[i] = i;
    for (size_t i = 0; i < t2.size(); ++i) t2[i] = i*i;
    ocv v2(t1.data(),t1.size());
    ocv v3(t2.data(),t2.size());
    ocv v4(v2);
    ocv v5(v3);
    /* Using formula Sum(i^3) = (N^2*(N+1)^2)/4 */
    double poly_sum = std::pow(t1.size()-1,2) * std::pow((t1.size()),2) / 4;
    result &= v2.dot(v3) == poly_sum; //First in memory; second in memory
    result &= v2.dot(v5) == poly_sum; //First in memory; second on disk
    result &= v4.dot(v3) == poly_sum; //First on disk; second in memory
    result &= v4.dot(v5) == poly_sum; //First on disk; second on disk
    ASSERT_TRUE(result);
}

TEST(OutOfCoreArray, dot_map) { //Also tests "==" operator
    using ocv = molpro::linalg::OutOfCoreArray<double, 100>;
    std::vector<double> v(10001);
    for (size_t i = 0; i < v.size(); ++i) v[i] = i;
    ocv v0(v.data(),v.size());
    ocv v1(v0);
    std::map<size_t,double> v_map;
    for (size_t i = 0; i < v.size(); ++i) { if (i%2 != 0) { v_map[i] = v[i]; } }
    bool result = true;
    result &= v0.dot(v_map) == v_map.size() * (2 * v_map.size() - 1) * (2 * v_map.size() + 1) / 3; //Vector in memory
    result &= v1.dot(v_map) == v_map.size() * (2 * v_map.size() - 1) * (2 * v_map.size() + 1) / 3; //Vector on disk
    ASSERT_TRUE(result);
}

TEST(OutOfCoreArray, axpy) {
    using ocv = molpro::linalg::OutOfCoreArray<double, 100>;
    std::vector<double> v(10001);
    for (size_t i = 0; i < v.size(); ++i) v[i] = 2 * i + 1;
    ocv v0(v.data(),v.size());
    ocv v1(v0);
    bool result = true;
    v0.axpy(2, v0); // Both point to ext. buffer
    v0.axpy(-3, v1); // One on disk
    result &= v0.dot(v0) < 1e-20;
    ASSERT_TRUE(result);
}

TEST(OutOfCoreArray, axpy_map) {
    using ocv = molpro::linalg::OutOfCoreArray<double, 100>;
    std::vector<double> v(10001);
    for (size_t i = 0; i < v.size(); ++i) v[i] = 2 * i + 1;
    ocv v0(v.data(),v.size());
    ocv v1(v0);
    std::map<size_t,double> v_map;
    for (size_t i = 0; i<50; i++) {
        size_t j = rand()%10001;
        v_map[j] = v1[j];
    }
    bool result = true;
    v0.axpy(2, v_map);
    v0.axpy(-2, v_map);
    v0.axpy(-1, v1);
    result &= v0.dot(v1) < 1e-20;
    ASSERT_TRUE(result);
}

TEST(OutOfCoreArray, scal) {
    using ocv = molpro::linalg::OutOfCoreArray<double, 100>;
    std::vector<double> v(10001);
    for (size_t i = 0; i < v.size(); ++i) v[i] = 2 * i + 1;
    ocv v0(v.data(),v.size());
    ocv v1(v0);
    bool result = true;
    v0.scal(3);
    v0.axpy(-3, v1);
    result &= v0.dot(v0) < 1e-20;
    ASSERT_TRUE(result);
}

TEST(OutOfCoreArray, zero) {
    using ocv = molpro::linalg::OutOfCoreArray<double, 100>;
    std::vector<double> v(10001);
    for (size_t i = 0; i < v.size(); ++i) v[i] = 2 * i + 1;
    ocv v0(v.data(),v.size());
    ocv v1(v0);
    bool result = true;
    v0.scal(0);
    result &= v0.dot(v0) < 1e-20;
    ASSERT_TRUE(result);
}

TEST(OutOfCoreArray, select) {
    using ocv = molpro::linalg::OutOfCoreArray<double, 100>;
    std::vector<double> v(10001);
    for (size_t i = 0; i < v.size(); ++i) v[i] = 2 * i + 1;
    ocv v0(v.data(), v.size());
    std::vector<double> t(10001);
    for (size_t i = 0; i < t.size(); ++i) t[i] = 1;
    ocv v1(t.data(), t.size());
    //ocv v1(v0);
    bool result = true;
    std::vector<size_t> indices;
    std::vector<double> values;
    std::tie(indices,values) = v0.select(v1);
    for (size_t i = 0; i < indices.size(); i++) {
        result &= v[indices[i]] == values[i]; // Check indices and values consistent
    }
    double max = *std::max_element(v.begin(),v.end());
    result &= max == v.back(); // Check largest value first
    ASSERT_TRUE(result);
}
