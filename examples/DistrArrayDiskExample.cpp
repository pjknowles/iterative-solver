#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/mpi.h>
#include <iostream>
int main(int argc, char* argv[]) {
  std::vector<double> v(10);
  for (size_t i = 0; i < v.size(); ++i)
    v[i] = i;
  molpro::mpi::init();
  molpro::linalg::array::Span<double> sv{v.data(), v.size()};
  auto [lo,hi]=
      molpro::linalg::array::util::make_distribution_spread_remainder<molpro::linalg::array::DistrArray::index_type>(
          v.size(), mpisize_global())
          .range(mpirank_global());
  molpro::linalg::array::DistrArraySpan dasv{v.size(), molpro::linalg::array::Span<double>(&sv[lo], hi - lo)};
  std::cout << dasv.dot(dasv)<<std::endl;
  molpro::linalg::array::DistrArrayFile dafv{v.size(),molpro::mpi::comm_global(),"/tmp"};
//  daf.put(0,v.size(),v.data());
  dafv.copy(dasv);
//  std::cout <<"f.f "<< dafv.dot(dafv)<<std::endl;
  std::cout <<"f.s "<< dafv.dot(dasv)<<std::endl;
//  std::cout <<"s.f "<< dasv.dot(dafv)<<std::endl;
  molpro::mpi::finalize();
  return 0;
}