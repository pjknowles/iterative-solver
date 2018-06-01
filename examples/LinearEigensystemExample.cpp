#include "IterativeSolver.h"
// Find lowest eigensolutions of M(i,j) = alpha*(i+1)*delta(i,j) + i + j
// Storage of vectors in-memory via class pv
using scalar = double;
constexpr size_t n = 200; // dimension of problem
constexpr scalar alpha = 100; // separation of diagonal elements

class pv : public LinearAlgebra::vector<scalar> {
public:
 std::vector<scalar> buffer;

 pv(size_t length = 0, int option = 0) : vector<scalar>(), buffer(length) {}

 pv(const pv *source, int option = 0) : vector<scalar>() { *this = *source; }

 void axpy(scalar a, const LinearAlgebra::vector<scalar> &other) override {
  for (size_t i = 0; i < buffer.size(); i++) buffer[i] += a * (dynamic_cast <const pv &> (other)).buffer[i];
 }

 void axpy(scalar a, const std::map<size_t, scalar> &other) override {
  for (const auto &o: other)
   (*this)[o.first] += a * o.second;
 }

 void scal(scalar a) override {
  for (size_t i = 0; i < buffer.size(); i++) buffer[i] *= a;
 }

 void zero() override {
  for (size_t i = 0; i < buffer.size(); i++) buffer[i] = 0;
 }

 scalar dot(const LinearAlgebra::vector<scalar> &other) const override {
  scalar result = 0;
  for (size_t i = 0; i < buffer.size(); i++) result += buffer[i] * (dynamic_cast <const pv &> (other)).buffer[i];
  return result;
 }

 scalar dot(const std::map<size_t, scalar> &other) const override {
  scalar result = 0;
  for (const auto &o: other)
   result += o.second * (*this)[o.first];
  return result;
 }

 std::tuple<std::vector<size_t>, std::vector<scalar> >
 select(const vector <scalar> &measure, const size_t maximumNumber = 1000, const scalar threshold = 0) const override {
  return std::make_tuple(std::vector<size_t>(0), std::vector<scalar>(0)); // null implementation
 };

 pv *clone(int option = 0) const override { return new pv(*this); }

 const scalar &operator[](const size_t i) const override { return buffer[i]; }

 scalar &operator[](size_t i) override { return buffer[i]; }
};

void action(const LinearAlgebra::vectorSet<scalar> &psx, LinearAlgebra::vectorSet<scalar> &outputs) {
 for (size_t k = 0; k < psx.size(); k++) {
  for (size_t i = 0; i < n; i++) {
   (*outputs[k])[i] = (alpha * (i + 1)) * (*psx[k])[i];
   for (size_t j = 0; j < n; j++)
    (*outputs[k])[i] += (i + j) * (*psx[k])[j];
  }
 }
}

void update(LinearAlgebra::vectorSet<scalar> &psc, const LinearAlgebra::vectorSet<scalar> &psg,
            std::vector<scalar> shift = std::vector<scalar>()) {
 for (size_t k = 0; k < psc.size(); k++)
  for (size_t i = 0; i < n; i++)
   (*psc[k])[i] = -(*psg[k])[i] / (shift[k] + 2 * i + alpha * (i + 1));
}

int main(int argc, char *argv[]) {
 LinearAlgebra::LinearEigensystem<scalar> solver;
 solver.m_verbosity = 1;
 solver.m_roots = 4;
 solver.m_thresh = 1e-6;
 LinearAlgebra::vectorSet<scalar> g;
 LinearAlgebra::vectorSet<scalar> x;
 for (int root = 0; root < solver.m_roots; root++) {
  x.push_back(std::make_shared<pv>(n));
  g.push_back(std::make_shared<pv>(n));
  x.back()->zero(); (*x.back())[root] = 1; // initial guess
 }
 std::vector<std::vector<scalar> > Pcoeff(solver.m_roots);
 for (auto iter = 0; iter < 1000; iter++) {
  action(x, g);
  solver.addVector(x, g, Pcoeff);
  update(x, g, solver.eigenvalues());
  if (solver.endIteration(x, g)) break;
 }
 std::cout << "Error={ ";
 for (int root = 0; root < solver.m_roots; root++)
  std::cout << solver.errors()[root] << " ";
 std::cout << "} after " << solver.iterations() << " iterations" << std::endl;
 for (int root = 0; root < solver.m_roots; root++) {
  std::cout << "Eigenvector:";
  for (size_t k = 0; k < n; k++) std::cout << " " << (*x[root])[k];
  std::cout << std::endl;
 }
}
