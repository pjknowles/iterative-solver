#include "IterativeSolver.h"
#include "PagedVector.h"

using namespace LinearAlgebra;

//  typedef SimpleParameterVector pv;
using scalar = double;
typedef PagedVector<scalar> pv;

static double alpha;
static double anharmonicity;
static double n;

void anharmonic_residual(const LinearAlgebra::vectorSet<scalar> &psx, LinearAlgebra::vectorSet<scalar> &outputs) {
  std::vector<scalar> psxk(n);
  std::vector<scalar> output(n);
  psx.front()->get(&(psxk[0]), n, 0);
  for (size_t i = 0; i < n; i++) {
    output[i] = (alpha * (i + 1) + anharmonicity * psxk[i]) * psxk[i];
    for (size_t j = 0; j < n; j++)
      output[i] += (i + j) * psxk[j];
  }
  outputs.front()->put(&output[0], n, 0);
}

void update(LinearAlgebra::vectorSet<scalar> &psc, const LinearAlgebra::vectorSet<scalar> &psg) {
  std::vector<scalar> psck(n);
  std::vector<scalar> psgk(n);
  psg.front()->get(&psgk[0], n, 0);
  psc.front()->get(&psck[0], n, 0);
  for (size_t i = 0; i < n; i++)
    psck[i] -= psgk[i] / (alpha * (i + 1));
  psc.front()->put(&psck[0], n, 0);
}

int main(int argc, char *argv[]) {
  alpha = 1;
  n = 100;
  anharmonicity = .5;
  DIIS<scalar> solver;
  solver.m_verbosity = 1;
  solver.m_maxIterations = 100;
  LinearAlgebra::vectorSet<scalar> g;
  LinearAlgebra::vectorSet<scalar> x;
  g.push_back(std::make_shared<pv>(n));
  x.push_back(std::make_shared<pv>(n));
  x.back()->zero();
  double one = 1;
  x.back()->put(&one, 1, 0);  // initial guess
  for (auto iter = 0; iter < solver.m_maxIterations; ++iter) {
    anharmonic_residual(x, g);
//     std::cout << "x: "<<x.front()<<std::endl;
//     std::cout << "g: "<<g.front()<<std::endl;
    solver.addVector(x, g);
    update(x, g);
    if (solver.endIteration(x, g)) break;
  }
//  if (not solver.solve(g,x)) std::cout << "Failure"<<std::endl;
  std::cout << "Distance of solution from origin: " << std::sqrt(x[0]->dot(*x[0])) << std::endl;
  std::cout << "Error=" << solver.errors().front() << " after " << solver.iterations() << " iterations" << std::endl;
}
