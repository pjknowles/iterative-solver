#include "IterativeSolver.h"
#include <cmath>
// Find lowest eigensolutions of M(i,j) = alpha*(i+1)*delta(i,j) + i + j
// Storage of vectors in-memory via class pv
using scalar = double;
constexpr size_t n = 300; // dimension of problem
constexpr scalar alpha = 100; // separation of diagonal elements

class pv {
 public:
  std::vector<scalar> buffer;
  using scalar_type = scalar;
  using element_type = scalar;

  explicit pv(size_t length = 0, int option = 0) : buffer(length) {}

  explicit pv(const pv& source, int option = 0) { *this = source; }

  void axpy(scalar a, const pv& other) {
    for (size_t i = 0; i < buffer.size(); i++) buffer[i] += a * other.buffer[i];
  }

  void axpy(scalar a, const std::map<size_t, scalar>& other) {
    for (const auto& o: other)
      (*this)[o.first] += a * o.second;
  }

  void scal(scalar a) {
    for (auto& b: buffer) b *= a;
  }

  void zero() {
    for (auto& b: buffer) b = 0;
  }

  scalar dot(const pv& other) const {
    scalar result = 0;
    for (size_t i = 0; i < buffer.size(); i++) result += buffer[i] * other.buffer[i];
    return result;
  }

  scalar dot(const std::map<size_t, scalar>& other) const {
    scalar result = 0;
    for (const auto& o: other)
      result += o.second * (*this)[o.first];
    return result;
  }

  std::tuple<std::vector<size_t>, std::vector<scalar> >
  select(const pv& measure, const size_t maximumNumber = 1000, const scalar threshold = 0) const {
    return std::make_tuple(std::vector<size_t>(0), std::vector<scalar>(0)); // null implementation
  };

  const scalar& operator[](const size_t i) const { return buffer[i]; }

  scalar& operator[](size_t i) { return buffer[i]; }
};

scalar matrix(const size_t i, const size_t j) {
  return (i == j ? alpha * (i + 1) : 0) + i + j;
}

void action(const std::vector<pv>& psx, std::vector<pv>& outputs) {
  for (size_t k = 0; k < psx.size(); k++) {
    for (size_t i = 0; i < n; i++) {
      outputs[k][i] = 0;
      for (size_t j = 0; j < n; j++)
        outputs[k][i] += matrix(i, j) * psx[k][j];
    }
  }
}

void update(std::vector<pv>& psc, const std::vector<pv>& psg,
            std::vector<scalar> shift = std::vector<scalar>()) {
  for (size_t k = 0; k < psc.size(); k++)
    for (size_t i = 0; i < n; i++)
      psc[k][i] -= psg[k][i] / (shift[k] + 2 * i + alpha * (i + 1));
}

int main(int argc, char* argv[]) {
  LinearAlgebra::LinearEigensystem<pv> solver;
  solver.m_verbosity = 1;
  solver.m_roots = 4;
  solver.m_thresh = 1e-6;
  std::vector<pv> g;
  std::vector<pv> x;
  std::vector<bool> active;
  for (size_t root = 0; root < solver.m_roots; root++) {
    active.push_back(true);
    x.emplace_back(n);
    g.emplace_back(n);
    x.back().zero();
    x.back()[root] = 1; // initial guess
  }
  std::vector<std::vector<scalar> > Pcoeff(solver.m_roots);
  for (auto iter = 0; iter < 1000; iter++) {
    action(x, g);
    solver.addVector(x, g, active, Pcoeff);
    update(x, g, solver.eigenvalues());
    if (solver.endIteration(x, g, active)) break;
  }
  std::cout << "Error={ ";
  for (const auto& e : solver.errors()) std::cout << e << " ";
  std::cout << "} after " << solver.iterations() << " iterations" << std::endl;
  for (size_t root = 0; root < solver.m_roots; root++) {
    std::cout << "Eigenvector: (norm=" << std::sqrt(x[root].dot(x[root])) << "): ";
    for (size_t k = 0; k < n; k++) std::cout << " " << (x[root])[k];
    std::cout << std::endl;
  }
}
