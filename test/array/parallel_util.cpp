#include "parallel_util.h"

namespace molpro {
namespace linalg {
namespace test {

std::tuple<std::vector<molpro::linalg::array::DistrArraySpan>, std::vector<molpro::linalg::array::DistrArraySpan>, std::vector<molpro::linalg::array::DistrArrayFile>> get_contiguous(size_t n, size_t dim) {

  //size_t n = 250;
  //size_t dim = 250;
  //std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim)),
  //    vz(n, std::vector<double>(dim));

  double* vx_mem = new double[n*dim];
  double* vy_mem = new double[n*dim];
  double* vz_mem = new double[n*dim];

//  std::vector<Span<double>> vx, vy, vz;
std::vector<molpro::linalg::array::span::Span<double>> vz;
//  vx.reserve(n);
//  vy.reserve(n);
  vz.reserve(n);
  for (size_t i=0; i<n; i++){
//    vx.emplace_back( Span(vx_mem+(i*n), dim) );
//    vy.emplace_back( Span(vy_mem+(i*n), dim) );
    vz.emplace_back( molpro::linalg::array::span::Span(vz_mem+(i*n), dim) );
  }

  std::vector<molpro::linalg::array::DistrArraySpan> cx;
  std::vector<molpro::linalg::array::DistrArraySpan> cy;
  std::vector<molpro::linalg::array::DistrArrayFile> cz;
  cx.reserve(n);
  cy.reserve(n);
  cz.reserve(n);
  
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(&vx_mem[i*dim], &vx_mem[(i*dim)+dim], i*100 + 0.5);
    std::iota(&vy_mem[i*dim], &vy_mem[(i*dim)+dim], 10000+ i*100 + 0.5);
    std::iota(&vz_mem[i*dim], &vz_mem[(i*dim)+dim], 20000+ i*100 + 0.5);
    cz.emplace_back(dim);
    auto crange = molpro::linalg::array::util::make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);

    auto clength = crange.second - crange.first;
    cx.emplace_back(dim, molpro::linalg::array::span::Span<double>(&vx_mem[i*dim + crange.first], clength), comm_global());
    cy.emplace_back(dim, molpro::linalg::array::span::Span<double>(&vy_mem[i*dim + crange.first], clength), comm_global());
//    std::cout << "i=" << i << " placing data in cz.back() range " << crange.first << ":" << crange.second << std::endl;
//    for (int j = 0; j < crange.second - crange.first; ++j)
//      std::cout << " " << vz_mem[i*dim+ crange.first + j];
//    std::cout << std::endl;
    cz.back().put(crange.first, crange.second, &vz_mem[i*dim + crange.first]);
  }

  return {cx, cy, cz};

}

} // namespace test
} // namespace linalg
} // namespace molpro