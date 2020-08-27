#include "molpro/linalg/IterativeSolver.h"
#include "molpro/linalg/SimpleArray.h"
#define USE_ARRAY 1
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif

#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/ArrayHandlerIterableSparse.h>
#include <molpro/linalg/array/ArrayHandlerSparse.h>

#include <algorithm>
#ifdef USE_ARRAY
#include <array>
#else
#include <vector>
#endif
#include <numeric>

using molpro::linalg::array::ArrayHandlerIterable;
using molpro::linalg::array::ArrayHandlerIterableSparse;
using molpro::linalg::array::ArrayHandlerSparse;
using molpro::linalg::iterativesolver::ArrayHandlers;

constexpr bool check = true;
constexpr bool print = false;

// For M(i,j) = alpha*(i+1)*delta(i,j) + i + j, b(i,n)=n+i
// solve M x = b
// Storage of vectors distributed and out of memory via SimpleArray class
using scalar = double;
// using pv = molpro::linalg::SimpleArray<scalar>;
constexpr size_t n = 300; // dimension of problem
#ifdef USE_ARRAY
using pv = std::array<scalar, n>;
#else
using pv = std::vector<scalar>;
#endif
using vectorSet = std::vector<pv>;
constexpr scalar alpha = 300; // separation of diagonal elements
// TODO nP>0
constexpr size_t nP = 0;    // number in initial P-space
constexpr size_t nRoot = 2; // number of equations

scalar matrix(const size_t i, const size_t j) { return (i == j ? alpha * (i + 1) : 0) + i + j; }

void action(const vectorSet& psx, vectorSet& outputs, size_t nroot) {
  for (size_t k = 0; k < nroot; k++) {
    for (size_t i = 0; i < n; i++) {
      outputs[k][i] = 0;
      for (size_t j = 0; j < n; j++)
        outputs[k][i] += matrix(i, j) * psx[k][j];
    }
  }
}

void actionP(const std::vector<std::map<size_t, scalar>> pspace, const std::vector<std::vector<scalar>>& psx,
             vectorSet& outputs, size_t nroot) {
  const size_t nP = pspace.size();
  for (size_t k = 0; k < nroot; k++) {
    for (size_t p = 0; p < nP; p++) {
      for (const auto& pc : pspace[p]) {
        for (size_t j = 0; j < n; j++)
          outputs[k][j] += matrix(pc.first, j) * pc.second * psx[k][p];
      }
    }
  }
}

void update(vectorSet& psc, const vectorSet& psg, size_t nroot) {
  for (size_t k = 0; k < nroot; k++) {
    for (size_t i = 0; i < n; i++)
      psc[k][i] -= psg[k][i] / (2 * i + alpha * (i + 1));
  }
}

int main(int argc, char* argv[]) {
  int rank = 0;
#ifdef HAVE_MPI_H
  MPI_Init(&argc, &argv);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    std::cout << "MPI size=" << size << std::endl;
#endif
  vectorSet g(nRoot);
  vectorSet b(nRoot);
  vectorSet x(nRoot);
  for (size_t root = 0; root < nRoot; root++) {
#ifndef USE_ARRAY
    g[root].resize(n);
    x[root].resize(n);
    b[root].resize(n);
#endif
    for (size_t i = 0; i < n; i++)
      b[root][i] = 1 / (1.0 + root + i);
  }
  std::vector<double> augmented_hessian_factors = {0, .001, .01, .1, 1};
  for (const auto& augmented_hessian_factor : augmented_hessian_factors) {
    std::cout << "Augmented hessian factor = " << augmented_hessian_factor << std::endl;
    auto handlers = std::make_shared<ArrayHandlers<pv, pv, std::map<size_t, double>>>();
    auto solver = molpro::linalg::LinearEquations<pv>(b, handlers, augmented_hessian_factor);
    solver.m_verbosity = 1;
    solver.m_roots = nRoot;
    solver.m_thresh = 1e-6;
    size_t p = 0;
    std::vector<scalar> PP;
    std::vector<std::map<size_t, scalar>> pspace(nP);
    for (auto& pc : pspace) {
      pc.clear();
      pc[p] = 1;
      for (size_t q = 0; q < nP; q++)
        PP.push_back(matrix(p, q));
      ++p;
    }
    std::vector<std::vector<scalar>> Pcoeff(solver.m_roots);
    for (size_t i = 0; i < solver.m_roots; ++i)
      Pcoeff[i].resize(nP);
    auto nwork = solver.addP(pspace, PP.data(), x, g, Pcoeff);
    for (auto iter = 0; iter < 100; iter++) {
      actionP(pspace, Pcoeff, g, nwork);
      update(x, g, nwork);
      if (rank == 0)
        solver.report();
      action(x, g, nwork);
      nwork = solver.addVector(x, g, Pcoeff);
      if (nwork == 0)
        break;
    }

    std::cout << "Error={ ";
    for (const auto& e : solver.errors())
      std::cout << e << " ";
    std::cout << "} after " << solver.statistics().iterations << " iterations" << std::endl;
    if (check or print) {
      std::vector<int> roots(solver.m_roots);
      std::iota(roots.begin(), roots.end(), 0);
      solver.solution(roots, x, g);
      if (print) {
        for (size_t root = 0; root < solver.m_roots; root++) {
          std::cout << "Solution:";
          for (size_t k = 0; k < n; k++)
            std::cout << " " << x[root][k];
          std::cout << std::endl;
        }
      }
      if (check) {
        action(x, g, solver.m_roots);
        for (size_t root = 0; root < solver.m_roots; root++) {
          std::transform(b[root].begin(), b[root].end(), g[root].begin(), g[root].begin(),
                         [](const scalar& a, const scalar& b) { return b - a; });
          auto err = std::inner_product(g[root].begin(), g[root].end(), g[root].begin(), 0);
          std::cout << "residual norm " << std::sqrt(err) << std::endl;
        }
      }
    }
  }
#ifdef HAVE_MPI_H
  MPI_Finalize();
#endif
}
