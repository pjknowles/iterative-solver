#include "test.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "molpro/linalg/IterativeSolver.h"
#include "molpro/linalg/PagedArray.h"
#include "molpro/linalg/SimpleArray.h"
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <limits>
#include <molpro/iostream.h>
#include <numeric>
#include <regex>
#include <vector>
// Find lowest eigensolutions of a matrix obtained from an external file
// Storage of vectors in-memory via class pv
using scalar = double;
size_t n;
std::vector<double> hmat;

class pv {
public:
  std::vector<scalar> buffer;
  using scalar_type = scalar;
  using value_type = scalar;

  explicit pv(size_t length = 0, int option = 0) : buffer(length) {}

  explicit pv(const pv& source, int option = 0) { *this = source; }

  void axpy(scalar a, const pv& other) {
    for (size_t i = 0; i < buffer.size(); i++)
      buffer[i] += a * other.buffer[i];
  }

  void axpy(scalar a, const std::map<size_t, scalar>& other) {
    for (const auto& o : other)
      (*this)[o.first] += a * o.second;
  }

  void scal(scalar a) {
    for (auto& b : buffer)
      b *= a;
  }

  void zero() {
    for (auto& b : buffer)
      b = 0;
  }

  scalar dot(const pv& other) const {
    scalar result = 0;
    for (size_t i = 0; i < buffer.size(); i++)
      result += buffer[i] * other.buffer[i];
    return result;
  }

  scalar dot(const std::map<size_t, scalar>& other) const {
    scalar result = 0;
    for (const auto& o : other)
      result += o.second * (*this)[o.first];
    return result;
  }

  std::tuple<std::vector<size_t>, std::vector<scalar>> select(const pv& measure, const size_t maximumNumber = 1000,
                                                              const scalar threshold = 0) const {
    return std::make_tuple(std::vector<size_t>(0), std::vector<scalar>(0)); // null implementation
  };

  const scalar& operator[](const size_t i) const { return buffer[i]; }

  scalar& operator[](size_t i) { return buffer[i]; }
};

scalar matrix(const size_t i, const size_t j) { return hmat[i * n + j]; }

void action(const std::vector<pv>& psx, std::vector<pv>& outputs) {
  for (size_t k = 0; k < psx.size(); k++) {
    for (size_t i = 0; i < n; i++) {
      outputs[k][i] = 0;
      for (size_t j = 0; j < n; j++)
        outputs[k][i] += matrix(i, j) * psx[k][j];
    }
  }
}

void update(std::vector<pv>& psc, const std::vector<pv>& psg, std::vector<scalar> shift = std::vector<scalar>()) {
  //  for (size_t k = 0; k < psc.size(); k++) {
  //   molpro::cout << "update root "<<k<<", shift="<<shift[k]<<": ";
  //    for (size_t i = 0; i < n; i++)
  //      molpro::cout <<" "<< -psg[k][i] / (1e-12-shift[k] + matrix(i,i));
  //    molpro::cout<<std::endl;
  //  }
  for (size_t k = 0; k < psc.size(); k++)
    for (size_t i = 0; i < n; i++)
      psc[k][i] -= psg[k][i] / (1e-12 - shift[k] + matrix(i, i));
}

TEST(IterativeSolver, file_eigen) {
  for (const auto& file : std::vector<std::string>{"bh"}) {
    std::vector<double> expected_eigenvalues;
    expected_eigenvalues.push_back(-25.127071042);
    expected_eigenvalues.push_back(-24.904605948);
    expected_eigenvalues.push_back(-24.858437717);
    expected_eigenvalues.push_back(-24.759563376);
    for (const auto& nroot : std::vector<int>{1, 2, 4}) {
      std::ifstream f(std::string{"./"} + file + ".hamiltonian");
      f >> n;
      molpro::cout << "\n\n*** " << file << ", " << nroot << " roots, problem dimension " << n << std::endl;
      hmat.resize(n * n);
      for (auto i = 0; i < n * n; i++)
        f >> hmat[i];
//      molpro::cout << "diagonal elements";
//      for (auto i = 0; i < n; i++)
//        molpro::cout << " " << matrix(i, i);
      std::vector<double> diagonals;
      for (auto i = 0; i < n; i++)
        diagonals.push_back(matrix(i, i));
      molpro::cout << std::endl;
      molpro::linalg::LinearEigensystem<pv> solver;
      solver.m_verbosity = 1;
      solver.m_roots = nroot;
      solver.m_thresh = 1e-9;
      std::vector<pv> g;
      std::vector<pv> x;
      for (size_t root = 0; root < solver.m_roots; root++) {
        x.emplace_back(n);
        g.emplace_back(n);
        x.back().zero();
        x.back()[std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin()] = 1; // initial guess
        *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
      }
      std::vector<std::vector<scalar>> Pcoeff(solver.m_roots);
      int nwork = solver.m_roots;
      for (auto iter = 0; iter < 100; iter++) {
        action(x, g);
        //        for (auto root = 0; root < nwork; root++) {
        //          std::cout << "before addVector() x:";
        //          for (auto i = 0; i < n; i++)
        //            std::cout << " " << x[root][i];
        //          std::cout << std::endl;
        //          std::cout << "before addVector() g:";
        //          for (auto i = 0; i < n; i++)
        //            std::cout << " " << g[root][i];
        //          std::cout << std::endl;
        //        }
        nwork = solver.addVector(x, g, Pcoeff);
        solver.report();
        if (nwork == 0)
          break;
        //                for (auto root = 0; root < nwork; root++) {
        //                  std::cout << "after addVector()"
        //                            << " eigenvalue=" << solver.working_set_eigenvalues()[root] << " error=" <<
        //                            solver.errors()[root]
        //                            << " x:";
        //                  for (auto i = 0; i < n; i++)
        //                    std::cout << " " << x[root][i];
        //                  std::cout << std::endl;
        //                  std::cout << "after addVector() g:";
        //                  for (auto i = 0; i < n; i++)
        //                    std::cout << " " << g[root][i];
        //                  std::cout << std::endl;
        //                }
        //        x.resize(nwork);
        //        g.resize(nwork);
        update(x, g, solver.working_set_eigenvalues());
        //            for (auto root = 0; root < nwork; root++) {
        //              std::cout << "after update() x:";
        //              for (auto i = 0; i < n; i++)
        //                std::cout << " " << x[root][i];
        //              std::cout << std::endl;
        //              std::cout << "after update() g:";
        //              for (auto i = 0; i < n; i++)
        //                std::cout << " " << g[root][i];
        //              std::cout << std::endl;
        //            }
      }
      std::cout << "Error={ ";
      for (const auto& e : solver.errors())
        std::cout << e << " ";
      std::cout << "} after " << solver.iterations() << " iterations" << std::endl;
      for (size_t root = 0; root < solver.m_roots; root++) {
        std::cout << "Eigenvalue " << std::fixed << std::setprecision(9) << solver.eigenvalues()[root] << std::endl;
        //        solver.solution(root, x.front(), g.front(), Pcoeff.front());
        //        std::cout << "Eigenvector: (norm=" << std::sqrt(x[0].dot(x[0])) << "): ";
        //        for (size_t k = 0; k < n; k++)
        //          std::cout << " " << (x[0])[k];
        //        std::cout << std::endl;
      }
      EXPECT_THAT(solver.errors(),
                  ::testing::Pointwise(::testing::DoubleNear(solver.m_thresh), std::vector<double>(nroot, double(0))));
      EXPECT_THAT(solver.eigenvalues(), ::testing::Pointwise(::testing::DoubleNear(2e-9),
                                                             std::vector<double>(expected_eigenvalues.data(), expected_eigenvalues.data() + solver.m_roots)));
    }
  }
}
