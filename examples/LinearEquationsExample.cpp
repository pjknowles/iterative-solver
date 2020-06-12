#include "molpro/linalg/IterativeSolver.h"
#include "molpro/linalg/array/PagedVector.h"
// For M(i,j) = alpha*(i+1)*delta(i,j) + i + j, b(i,n)=n+i
// solve M x = b
// Storage of vectors distributed and out of memory via PagedVector class
using scalar = double;
using pv = molpro::linalg::PagedVector<scalar>;
using vectorSet = std::vector<pv>;
constexpr size_t n = 300; // dimension of problem
constexpr scalar alpha = 300; // separation of diagonal elements
constexpr size_t nP = 10; // number in initial P-space
constexpr size_t nRoot = 2; // number of equations

scalar matrix(const size_t i, const size_t j) {
  return (i == j ? alpha * (i + 1) : 0) + i + j;
}

void action(const vectorSet& psx, vectorSet& outputs) {
  std::vector<scalar> psxk(n);
  std::vector<scalar> output(n);
  for (size_t k = 0; k < psx.size(); k++) {
    psx[k].get(&(psxk[0]), n, 0);
    for (size_t i = 0; i < n; i++) {
      output[i] = 0;
      for (size_t j = 0; j < n; j++)
        output[i] += matrix(i, j) * psxk[j];
    }
    outputs[k].put(&output[0], n, 0);
  }
}

void actionP(
    const std::vector<std::map<size_t, scalar> > pspace,
    const std::vector<std::vector<scalar> >& psx,
    vectorSet& outputs) {
  const size_t nP = pspace.size();
  for (size_t k = 0; k < psx.size(); k++) {
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

void update(vectorSet& psc, const vectorSet& psg) {
  std::vector<scalar> psck(n);
  std::vector<scalar> psgk(n);
  for (size_t k = 0; k < psc.size(); k++) {
    psg[k].get(&psgk[0], n, 0);
    psc[k].get(&psck[0], n, 0);
    for (size_t i = 0; i < n; i++)
      psck[i] -= psgk[i] / (2 * i + alpha * (i + 1));
    psc[k].put(&psck[0], n, 0);
  }
}

int main(int argc, char* argv[]) {
  vectorSet g; g.reserve(nRoot);
  vectorSet b; b.reserve(nRoot);
  vectorSet x; x.reserve(nRoot);
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
    molpro::linalg::LinearEquations<pv> solver(b, augmented_hessian_factor);
    solver.m_verbosity = 1;
    solver.m_roots = nRoot;
    solver.m_thresh = 1e-6;
    size_t p = 0;
    std::vector<scalar> PP;
    std::vector<std::map<size_t, scalar> > pspace(nP);
    for (auto& pc : pspace) {
      pc.clear();
      pc[p] = 1;
      for (size_t q = 0; q < nP; q++) PP.push_back(matrix(p, q));
      ++p;
    }
    std::vector<std::vector<scalar> > Pcoeff(solver.m_roots);
    for (size_t i = 0; i < solver.m_roots; ++i) Pcoeff[i].resize(nP);
    solver.addP(pspace, PP.data(), x, g, Pcoeff);
    for (auto iter = 0; iter < 100; iter++) {
      actionP(pspace, Pcoeff, g);
      update(x, g);
      if (solver.endIteration(x, g)) break;
      action(x, g);
      solver.addVector(x, g, Pcoeff);
    }

    std::cout << "Error={ ";
    for (const auto& e : solver.errors()) std::cout << e << " ";
    std::cout << "} after " << solver.iterations() << " iterations" << std::endl;
    action(x, g);
    for (size_t root = 0; root < solver.m_roots; root++) {
      g[root].axpy(-1, b[root]);
      std::vector<scalar> buf(n);
      x[root].get(&buf[0], n, 0);
      std::cout << "Solution:";
      for (size_t k = 0; k < n; k++) std::cout << " " << buf[k];
      std::cout << std::endl;
      auto err = g[root].dot(g[root]);
      std::cout << "residual norm " << std::sqrt(err) << std::endl;
    }
  }
}
