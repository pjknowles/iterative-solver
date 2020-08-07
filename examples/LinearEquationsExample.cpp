#include "molpro/linalg/IterativeSolver.h"
#include "molpro/linalg/SimpleArray.h"
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif

#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/ArrayHandlerIterableSparse.h>
#include <molpro/linalg/array/ArrayHandlerSparse.h>

using molpro::linalg::array::ArrayHandlerIterable;
using molpro::linalg::array::ArrayHandlerIterableSparse;
using molpro::linalg::array::ArrayHandlerSparse;
using molpro::linalg::iterativesolver::ArrayHandlers;

// For M(i,j) = alpha*(i+1)*delta(i,j) + i + j, b(i,n)=n+i
// solve M x = b
// Storage of vectors distributed and out of memory via SimpleArray class
using scalar = double;
using pv = molpro::linalg::SimpleArray<scalar>;
using vectorSet = std::vector<pv>;
constexpr size_t n = 300;     // dimension of problem
constexpr scalar alpha = 300; // separation of diagonal elements
// TODO nP>0
constexpr size_t nP = 0;    // number in initial P-space
constexpr size_t nRoot = 2; // number of equations

scalar matrix(const size_t i, const size_t j) { return (i == j ? alpha * (i + 1) : 0) + i + j; }

void action(const vectorSet& psx, vectorSet& outputs, size_t nroot) {
  std::vector<scalar> psxk(n);
  std::vector<scalar> output(n);
  for (size_t k = 0; k < nroot; k++) {
    psx[k].get(&(psxk[0]), n, 0);
    for (size_t i = 0; i < n; i++) {
      output[i] = 0;
      for (size_t j = 0; j < n; j++)
        output[i] += matrix(i, j) * psxk[j];
    }
    outputs[k].put(&output[0], n, 0);
  }
}

void actionP(const std::vector<std::map<size_t, scalar>> pspace, const std::vector<std::vector<scalar>>& psx,
             vectorSet& outputs, size_t nroot) {
  const size_t nP = pspace.size();
  for (size_t k = 0; k < nroot; k++) {
    std::vector<scalar> output(n);
    outputs[k].get(&output[0], n, 0);
    for (size_t p = 0; p < nP; p++) {
      for (const auto& pc : pspace[p]) {
        for (size_t j = 0; j < n; j++)
          output[j] += matrix(pc.first, j) * pc.second * psx[k][p];
      }
    }
    outputs[k].put(&output[0], n, 0);
  }
}

void update(vectorSet& psc, const vectorSet& psg, size_t nroot) {
  std::vector<scalar> psck(n);
  std::vector<scalar> psgk(n);
  for (size_t k = 0; k < nroot; k++) {
    psg[k].get(&psgk[0], n, 0);
    psc[k].get(&psck[0], n, 0);
    for (size_t i = 0; i < n; i++)
      psck[i] -= psgk[i] / (2 * i + alpha * (i + 1));
    psc[k].put(&psck[0], n, 0);
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
  vectorSet g;
  g.reserve(nRoot);
  vectorSet b;
  b.reserve(nRoot);
  vectorSet x;
  x.reserve(nRoot);
  for (size_t root = 0; root < nRoot; root++) {
    std::vector<scalar> bb(n);
    for (size_t i = 0; i < n; i++)
      bb[i] = 1 / (1.0 + root + i);
    b.emplace_back(n);
    b.back().put(bb.data(), bb.size(), 0);
    x.emplace_back(n);
    g.emplace_back(n);
  }
  std::vector<double> augmented_hessian_factors = {0, .001, .01, .1, 1};
  for (const auto& augmented_hessian_factor : augmented_hessian_factors) {
    std::cout << "Augmented hessian factor = " << augmented_hessian_factor << std::endl;
    auto rr = std::make_shared<ArrayHandlerIterable<pv>>();
    auto qq = std::make_shared<ArrayHandlerIterable<pv>>();
    auto pp = std::make_shared<ArrayHandlerSparse<std::map<size_t, double>>>();
    auto rq = std::make_shared<ArrayHandlerIterable<pv>>();
    auto rp = std::make_shared<ArrayHandlerIterableSparse<pv, std::map<size_t, double>>>();
    auto qr = std::make_shared<ArrayHandlerIterable<pv>>();
    auto qp = std::make_shared<ArrayHandlerIterableSparse<pv, std::map<size_t, double>>>();
    auto handlers = ArrayHandlers<pv, pv, std::map<size_t, double>>{rr, qq, pp, rq, rp, qr, qp};
    molpro::linalg::LinearEquations<pv> solver(b, handlers, augmented_hessian_factor);
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
    std::cout << "} after " << solver.iterations() << " iterations" << std::endl;
    if (false) {
      std::vector<int> roots(solver.m_roots);
      std::iota(roots.begin(), roots.end(), 0);
      solver.solution(roots, x, g);
      action(x, g, solver.m_roots);
      for (size_t root = 0; root < solver.m_roots; root++) {
        g[root].axpy(-1, b[root]);
        std::vector<scalar> buf(n);
        x[root].get(&buf[0], n, 0);
        std::cout << "Solution:";
        for (size_t k = 0; k < n; k++)
          std::cout << " " << buf[k];
        std::cout << std::endl;
        auto err = g[root].dot(g[root]);
        std::cout << "residual norm " << std::sqrt(err) << std::endl;
      }
    }
  }
}
