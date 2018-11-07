#include "IterativeSolver.h"
#include "PagedVector.h"


//  typedef SimpleParameterVector pv;
using scalar = double;
using pv = LinearAlgebra::PagedVector<scalar>;

static double alpha;
static double anharmonicity;
static double n;

void anharmonic_residual(const std::vector<pv>& psx, std::vector<pv>& outputs) {
  std::vector<scalar> psxk(n);
  std::vector<scalar> output(n);
  psx.front().get(psxk.data(), n, 0);
  for (size_t i = 0; i < n; i++) {
    output[i] = (alpha * (i + 1) + anharmonicity * psxk[i]) * psxk[i];
    for (size_t j = 0; j < n; j++)
      output[i] += (i + j) * psxk[j];
  }
  outputs.front().put(output.data(), n, 0);
}

void update(std::vector<pv>& psc, const std::vector<pv>& psg) {
  std::vector<scalar> psck(n);
  std::vector<scalar> psgk(n);
  psg.front().get(psgk.data(), n, 0);
  psc.front().get(psck.data(), n, 0);
  for (size_t i = 0; i < n; i++)
    psck[i] -= psgk[i] / (alpha * (i + 1));
  psc.front().put(psck.data(), n, 0);
}

int main(int argc, char* argv[]) {
  alpha = 1;
  n = 100;
  anharmonicity = .5;
  LinearAlgebra::DIIS<pv> solver;
  solver.m_verbosity = 1;
  solver.m_maxIterations = 100;
  std::vector<pv> g;
  std::vector<pv> x;
  std::vector<bool> active;
  g.emplace_back(n);
  x.emplace_back(n);
  x.back().scal(0);
  scalar one = 1;
  x.back().put(&one, 1, 0);  // initial guess
  active.push_back(true);
  for (size_t iter = 0; iter < solver.m_maxIterations; ++iter) {
    anharmonic_residual(x, g);
    solver.addVector(x, g, active);
    update(x, g);
    if (solver.endIteration(x, g, active)) break;
  }
  std::cout << "Distance of solution from origin: " << std::sqrt(x[0].dot(x[0])) << std::endl;
  std::cout << "Error=" << solver.errors().front() << " after " << solver.iterations() << " iterations" << std::endl;
}
