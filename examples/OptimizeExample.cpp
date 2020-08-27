#include "molpro/linalg/IterativeSolver.h"
#include "molpro/linalg/SimpleArray.h"
#include <regex>

#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/ArrayHandlerIterableSparse.h>
#include <molpro/linalg/array/ArrayHandlerSparse.h>

using molpro::linalg::array::ArrayHandlerIterable;
using molpro::linalg::array::ArrayHandlerIterableSparse;
using molpro::linalg::array::ArrayHandlerSparse;
using molpro::linalg::iterativesolver::ArrayHandlers;

//  typedef SimpleParameterVector pv;
using scalar = double;
using pv = molpro::linalg::SimpleArray<scalar>;

static double alpha;
static double anharmonicity;
static size_t n;

scalar anharmonic_residual(const pv& psx, pv& outputs) {
  std::vector<scalar> psxk(n);
  std::vector<scalar> output(n);
  psx.get(psxk.data(), n, 0);
  scalar value = 0;
  for (size_t i = 0; i < n; i++) {
    value += (alpha * (i + 1) / scalar(2) + anharmonicity * (psxk[i] - 1) / scalar(3)) * (psxk[i] - 1) * (psxk[i] - 1);
    output[i] = (alpha * (i + 1) + anharmonicity * (psxk[i] - 1)) * (psxk[i] - 1);
    for (size_t j = 0; j < n; j++)
      output[i] += (i + j) * (psxk[j] - 1);
  }
  outputs.put(output.data(), n, 0);
  return value;
}

void update(pv& psc, const pv& psg) {
  std::vector<scalar> psck(n);
  std::vector<scalar> psgk(n);
  psg.get(psgk.data(), n, 0);
  psc.get(psck.data(), n, 0);
  for (size_t i = 0; i < n; i++)
    psck[i] -= psgk[i] / (2 * i + alpha * (i + 1));
  psc.put(psck.data(), n, 0);
}

int main(int argc, char* argv[]) {
  alpha = 7;
  n = 10;
  anharmonicity = 0.7;
  for (const auto& method : std::vector<std::string>{"null", "L-BFGS"}) {
    std::cout << "optimize with " << method << std::endl;
    auto handlers = std::make_shared<ArrayHandlers<pv, pv, std::map<size_t, double>>>();
    molpro::linalg::Optimize<pv> solver(handlers, std::regex_replace(method, std::regex("-iterate"), ""));
    solver.m_maxIterations = 100;
    solver.m_verbosity = 1;
    pv g(n);
    pv x(n);
    pv hg(n);
    scalar one = 1;
    for (auto i = 0; i < n; i++)
      x.put(&one, 1, i);
    scalar zero = 0;
    x.put(&zero, 1, 0); // initial guess
    for (size_t iter = 1; iter <= solver.m_maxIterations; ++iter) {
      auto value = anharmonic_residual(x, g);
      //      xout << "iteration "<<iter<<" value="<<value<<"\n x: "<<x<<"\n g: "<<g<<std::endl;
      if (solver.addValue(x, value, g))
        update(x, g);
      if (solver.endIteration(x, g))
        break;
    }
    for (size_t i = 0; i < x.size(); i++) {
      scalar val;
      x.get(&val, 1, i);
      val -= 1;
      x.put(&val, 1, i);
    }
    std::cout << "Distance of solution from exact solution: " << std::sqrt(x.dot(x)) << std::endl;
    std::cout << "Error=" << solver.errors().front() << " after " << solver.statistics().iterations << " iterations"
              << std::endl;
  }
}
