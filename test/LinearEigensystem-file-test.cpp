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
#include <molpro/linalg/SimpleArray.h>
#include <numeric>
#include <regex>
#include <vector>

#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/ArrayHandlerIterableSparse.h>
#include <molpro/linalg/array/ArrayHandlerSparse.h>

using molpro::linalg::array::ArrayHandler;
using molpro::linalg::array::ArrayHandlerIterable;
using molpro::linalg::array::ArrayHandlerIterableSparse;
using molpro::linalg::array::ArrayHandlerSparse;
using molpro::linalg::itsolv::ArrayHandlers;
// Find lowest eigensolution of a matrix obtained from an external file
// Storage of vectors in-memory via class pv
using scalar = double;
size_t n;
std::vector<double> hmat;

using pv = molpro::linalg::SimpleArray<scalar>;

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
std::vector<double> residual(const std::vector<pv>& psx, std::vector<pv>& outputs) {

  std::vector<double> evals;
  for (size_t k = 0; k < psx.size(); k++) {
    double eval = 0;
    for (size_t i = 0; i < n; i++) {
      outputs[k][i] = 0;
      for (size_t j = 0; j < n; j++)
        outputs[k][i] += matrix(i, j) * psx[k][j];
      eval += outputs[k][i] * psx[k][i];
    }
    eval /= psx[k].dot(psx[k]);
    for (size_t i = 0; i < n; i++) {
      outputs[k][i] -= eval * psx[k][i];
    }
    evals.push_back(eval);
  }
  return evals;
}

void update(std::vector<pv>& psc, const std::vector<pv>& psg, std::vector<scalar> shift = std::vector<scalar>()) {
  //  for (size_t k = 0; k < psc.size(); k++) {
  //   molpro::cout << "update root "<<k<<", shift="<<shift[k]<<": ";
  //    for (size_t i = 0; i < n; i++)
  //      molpro::cout <<" "<< -psg[k][i] / (1e-12-shift[k] + matrix(i,i));
  //    molpro::cout<<std::endl;
  //  }
  for (size_t k = 0; k < psc.size() && k < shift.size(); k++)
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
    for (const auto& nroot : std::vector<int>{1, 2, 3}) {
      std::ifstream f(std::string{"./"} + file + ".hamiltonian");
      f >> n;
      hmat.resize(n * n);
      for (auto i = 0; i < n * n; i++)
        f >> hmat[i];
      //      molpro::cout << "diagonal elements";
      //      for (auto i = 0; i < n; i++)
      //        molpro::cout << " " << matrix(i, i);
      //      molpro::cout << std::endl;
      for (auto np = 0; np <= 16; np += 4) {
        molpro::cout << "\n\n*** " << file << ", " << nroot << " roots, problem dimension " << n
                     << ", pspace dimension " << np << std::endl;
        auto handlers = std::make_shared<ArrayHandlers<pv, pv, std::map<size_t, double>>>();
        auto solver = molpro::linalg::LinearEigensystem<pv>{handlers};
        solver.m_verbosity = 1;
        solver.m_roots = nroot;
        solver.m_thresh = 1e-9;
        std::vector<pv> g;
        std::vector<pv> x;
        std::vector<size_t> guess;
        {
          std::vector<double> diagonals;
          for (auto i = 0; i < n; i++)
            diagonals.push_back(matrix(i, i));
          for (size_t root = 0; root < solver.m_roots; root++) {
            x.emplace_back(n);
            g.emplace_back(n);
            x.back().scal(0);
            guess.push_back(std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin()); // initial guess
            *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
            x.back()[guess.back()] = 1; // initial guess
          }
        }
        //        std::cout << "guess: " << guess << std::endl;
        int nwork = solver.m_roots;

        std::vector<std::vector<scalar>> Pcoeff(solver.m_roots);
        std::vector<scalar> PP;
        std::vector<std::map<size_t, scalar>> pspace;

        for (auto iter = 0; iter < 32; iter++) {
          if (iter == 0 && np > 0) {
            std::vector<double> diagonals;
            for (auto i = 0; i < n; i++)
              diagonals.push_back(matrix(i, i));
            for (size_t p = 0; p < np; p++) {
              std::map<size_t, scalar> pp;
              pp.insert({std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin(), 1});
              pspace.push_back(pp);
              *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
            }
            //            std::cout << "P space: " << pspace << std::endl;
            PP.reserve(np * np);
            for (const auto& i : pspace)
              for (const auto& j : pspace)
                PP.push_back(matrix(i.begin()->first, j.begin()->first));
            //            std::cout << "PP: " << PP << std::endl;
            nwork = solver.addP(pspace, PP.data(), x, g, Pcoeff);
          } else {
            action(x, g);
            //            for (auto root = 0; root < nwork; root++) {
            //              std::cout << "before addVector() x:";
            //              for (auto i = 0; i < n; i++)
            //                std::cout << " " << x[root][i];
            //              std::cout << std::endl;
            //              std::cout << "before addVector() g:";
            //              for (auto i = 0; i < n; i++)
            //                std::cout << " " << g[root][i];
            //              std::cout << std::endl;
            //            }
            nwork = solver.addVector(x, g, Pcoeff);
          }
          if (nwork == 0)
            break;
          //          for (auto root = 0; root < nwork; root++) {
          //            std::cout << "after addVector()"
          //                      << " eigenvalue=" << solver.working_set_eigenvalues()[root] << " error=" <<
          //                      solver.errors()[root]
          //                      << " x:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << x[root][i];
          //            std::cout << std::endl;
          //            std::cout << "after addVector() g:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << g[root][i];
          //            std::cout << std::endl;
          //          }
          //        x.resize(nwork);
          //        g.resize(nwork);
          for (auto k = 0; k < solver.working_set().size(); k++) {
            for (auto p = 0; p < np; p++)
              for (const auto& pc : pspace[p])
                for (auto i = 0; i < n; i++)
                  g[k][i] += matrix(i, pc.first) * pc.second * Pcoeff[k][p];
            //            molpro::cout << "old error " << solver.m_errors[solver.working_set()[k]];
            solver.m_errors[solver.working_set()[k]] = std::sqrt(g[k].dot(g[k]));
            //            molpro::cout << "; new error " << solver.m_errors[solver.working_set()[k]] << std::endl;
          }
          //          for (auto root = 0; root < nwork; root++) {
          //            std::cout << "after p gradient g:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << g[root][i];
          //            std::cout << std::endl;
          //          }
          update(x, g, solver.working_set_eigenvalues());
          //          for (auto root = 0; root < nwork; root++) {
          //            std::cout << "after update() x:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << x[root][i];
          //            std::cout << std::endl;
          //            std::cout << "after update() g:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << g[root][i];
          //            std::cout << std::endl;
          //          }
          solver.report();
          if (*std::max_element(solver.errors().begin(), solver.errors().end()) < solver.m_thresh)
            break;
        }
        std::cout << "Error={ ";
        for (const auto& e : solver.errors())
          std::cout << e << " ";
        std::cout << "} after " << solver.statistics().iterations << " iterations" << std::endl;
        for (size_t root = 0; root < solver.m_roots; root++) {
          std::cout << "Eigenvalue " << std::fixed << std::setprecision(9) << solver.eigenvalues()[root] << std::endl;
          //        solver.solution(root, x.front(), g.front(), Pcoeff.front());
          //        std::cout << "Eigenvector: (norm=" << std::sqrt(x[0].dot(x[0])) << "): ";
          //        for (size_t k = 0; k < n; k++)
          //          std::cout << " " << (x[0])[k];
          //        std::cout << std::endl;
        }
        EXPECT_THAT(solver.errors(), ::testing::Pointwise(::testing::DoubleNear(2*solver.m_thresh),
                                                          std::vector<double>(nroot, double(0))));
        EXPECT_THAT(solver.eigenvalues(),
                    ::testing::Pointwise(::testing::DoubleNear(2e-9),
                                         std::vector<double>(expected_eigenvalues.data(),
                                                             expected_eigenvalues.data() + solver.m_roots)));
      }
    }
  }
}

TEST(IterativeSolver, file_optimize_eigenvalue) {
  for (const auto& file : std::vector<std::string>{"bh"}) {
    std::vector<double> expected_eigenvalues;
    expected_eigenvalues.push_back(-25.127071042);
    expected_eigenvalues.push_back(-24.904605948);
    expected_eigenvalues.push_back(-24.858437717);
    expected_eigenvalues.push_back(-24.759563376);
    {
      std::ifstream f(std::string{"./"} + file + ".hamiltonian");
      f >> n;
      hmat.resize(n * n);
      for (auto i = 0; i < n * n; i++)
        f >> hmat[i];
      //      molpro::cout << "diagonal elements";
      //      for (auto i = 0; i < n; i++)
      //        molpro::cout << " " << matrix(i, i);
      //      molpro::cout << std::endl;
      {
        molpro::cout << "\n\n*** " << file << ", BFGS for lowest eigensolution, problem dimension " << n << std::endl;
        auto rr = std::make_shared<ArrayHandlerIterable<pv>>();
        auto qq = std::make_shared<ArrayHandlerIterable<pv>>();
        auto pp = std::make_shared<ArrayHandlerSparse<std::map<size_t, double>>>();
        auto rq = std::make_shared<ArrayHandlerIterable<pv>>();
        auto rp = std::make_shared<ArrayHandlerIterableSparse<pv, std::map<size_t, double>>>();
        auto qr = std::make_shared<ArrayHandlerIterable<pv>>();
        auto qp = std::make_shared<ArrayHandlerIterableSparse<pv, std::map<size_t, double>>>();
        auto handlers = std::make_shared<ArrayHandlers<pv, pv, std::map<size_t, double>>>(rr, qq, pp, rq, rp, qr, qp);
        auto solver = molpro::linalg::Optimize<pv>{handlers};
        solver.m_verbosity = 1;
        solver.m_thresh = 1e-9;
        std::vector<pv> g;
        std::vector<pv> x;
        std::vector<size_t> guess;
        {
          std::vector<double> diagonals;
          for (auto i = 0; i < n; i++)
            diagonals.push_back(matrix(i, i));
          for (size_t root = 0; root < solver.m_roots; root++) {
            x.emplace_back(n);
            g.emplace_back(n);
            x.back().scal(0);
            guess.push_back(std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin()); // initial guess
            *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
            x.back()[guess.back()] = 1; // initial guess
          }
        }
        double e0 = matrix(guess.back(), guess.back());
        //        std::cout << "guess: " << guess << std::endl;
        int nwork = 1;

        for (auto iter = 0; iter < 50; iter++) {
          auto eval = residual(x, g).front();
          //            for (auto root = 0; root < nwork; root++) {
          //              std::cout << "before addVector() x:";
          //              for (auto i = 0; i < n; i++)
          //                std::cout << " " << x[root][i];
          //              std::cout << std::endl;
          //              std::cout << "before addVector() g:";
          //              for (auto i = 0; i < n; i++)
          //                std::cout << " " << g[root][i];
          //              std::cout << std::endl;
          //            }
          nwork = solver.addValue(x.front(), eval, g.front());
          if (nwork == 0)
            break;
          //          for (auto root = 0; root < nwork; root++) {
          //            std::cout << "after addVector()"
          //                      << " eigenvalue=" << solver.working_set_eigenvalues()[root] << " error=" <<
          //                      solver.errors()[root]
          //                      << " x:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << x[root][i];
          //            std::cout << std::endl;
          //            std::cout << "after addVector() g:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << g[root][i];
          //            std::cout << std::endl;
          //          }
          //        x.resize(nwork);
          //        g.resize(nwork);
          //          }
          update(x, g, {eval});
          //          for (auto root = 0; root < nwork; root++) {
          //            std::cout << "after update() x:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << x[root][i];
          //            std::cout << std::endl;
          //            std::cout << "after update() g:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << g[root][i];
          //            std::cout << std::endl;
          //          }
          solver.endIteration(x, g);
          if (*std::max_element(solver.errors().begin(), solver.errors().end()) < solver.m_thresh)
            break;
        }
        std::cout << "Error={ ";
        for (const auto& e : solver.errors())
          std::cout << e << " ";
        std::cout << "} after " << solver.statistics().iterations << " iterations" << std::endl;
        auto evals = residual(x, g);
        for (size_t root = 0; root < solver.m_roots; root++) {
          std::cout << "Eigenvalue " << std::fixed << std::setprecision(9) << evals[root] << std::endl;
          //        solver.solution(root, x.front(), g.front(), Pcoeff.front());
          //        std::cout << "Eigenvector: (norm=" << std::sqrt(x[0].dot(x[0])) << "): ";
          //        for (size_t k = 0; k < n; k++)
          //          std::cout << " " << (x[0])[k];
          //        std::cout << std::endl;
        }
        EXPECT_THAT(solver.errors(),
                    ::testing::Pointwise(::testing::DoubleNear(solver.m_thresh), std::vector<double>(1, double(0))));
        EXPECT_THAT(evals, ::testing::Pointwise(::testing::DoubleNear(2e-9),
                                                std::vector<double>(expected_eigenvalues.data(),
                                                                    expected_eigenvalues.data() + solver.m_roots)));
      }
    }
  }
}

TEST(IterativeSolver, file_diis_eigenvalue) {
  for (const auto& file : std::vector<std::string>{"bh"}) {
    std::vector<double> expected_eigenvalues;
    expected_eigenvalues.push_back(-25.127071042);
    expected_eigenvalues.push_back(-24.904605948);
    expected_eigenvalues.push_back(-24.858437717);
    expected_eigenvalues.push_back(-24.759563376);
    {
      std::ifstream f(std::string{"./"} + file + ".hamiltonian");
      f >> n;
      hmat.resize(n * n);
      for (auto i = 0; i < n * n; i++)
        f >> hmat[i];
      //      molpro::cout << "diagonal elements";
      //      for (auto i = 0; i < n; i++)
      //        molpro::cout << " " << matrix(i, i);
      //      molpro::cout << std::endl;
      {
        molpro::cout << "\n\n*** " << file << ", DIIS for lowest eigensolution, problem dimension " << n << std::endl;
        auto rr = std::make_shared<ArrayHandlerIterable<pv>>();
        auto qq = std::make_shared<ArrayHandlerIterable<pv>>();
        auto pp = std::make_shared<ArrayHandlerSparse<std::map<size_t, double>>>();
        auto rq = std::make_shared<ArrayHandlerIterable<pv>>();
        auto rp = std::make_shared<ArrayHandlerIterableSparse<pv, std::map<size_t, double>>>();
        auto qr = std::make_shared<ArrayHandlerIterable<pv>>();
        auto qp = std::make_shared<ArrayHandlerIterableSparse<pv, std::map<size_t, double>>>();
        auto handlers = std::make_shared<ArrayHandlers<pv, pv, std::map<size_t, double>>>(rr, qq, pp, rq, rp, qr, qp);
        auto solver = molpro::linalg::DIIS<pv>{handlers};
        solver.m_verbosity = 1;
        solver.m_thresh = 1e-9;
        std::vector<pv> g;
        std::vector<pv> x;
        std::vector<size_t> guess;
        {
          std::vector<double> diagonals;
          for (auto i = 0; i < n; i++)
            diagonals.push_back(matrix(i, i));
          for (size_t root = 0; root < solver.m_roots; root++) {
            x.emplace_back(n);
            g.emplace_back(n);
            x.back().scal(0);
            guess.push_back(std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin()); // initial guess
            *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
            x.back()[guess.back()] = 1; // initial guess
          }
        }
        double e0 = matrix(guess.back(), guess.back());
        //        std::cout << "guess: " << guess << std::endl;
        int nwork = 1;

        for (auto iter = 0; iter < 50; iter++) {
          auto eval = residual(x, g).front();
          //            for (auto root = 0; root < nwork; root++) {
          //              std::cout << "before addVector() x:";
          //              for (auto i = 0; i < n; i++)
          //                std::cout << " " << x[root][i];
          //              std::cout << std::endl;
          //              std::cout << "before addVector() g:";
          //              for (auto i = 0; i < n; i++)
          //                std::cout << " " << g[root][i];
          //              std::cout << std::endl;
          //            }
          nwork = solver.addVector(x.front(), g.front());
          if (nwork == 0)
            break;
          //          for (auto root = 0; root < nwork; root++) {
          //            std::cout << "after addVector()"
          //                      << " eigenvalue=" << solver.working_set_eigenvalues()[root] << " error=" <<
          //                      solver.errors()[root]
          //                      << " x:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << x[root][i];
          //            std::cout << std::endl;
          //            std::cout << "after addVector() g:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << g[root][i];
          //            std::cout << std::endl;
          //          }
          //        x.resize(nwork);
          //        g.resize(nwork);
          //          }
          update(x, g, {eval});
          //          for (auto root = 0; root < nwork; root++) {
          //            std::cout << "after update() x:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << x[root][i];
          //            std::cout << std::endl;
          //            std::cout << "after update() g:";
          //            for (auto i = 0; i < n; i++)
          //              std::cout << " " << g[root][i];
          //            std::cout << std::endl;
          //          }
          solver.endIteration(x, g);
          if (*std::max_element(solver.errors().begin(), solver.errors().end()) < solver.m_thresh)
            break;
        }
        std::cout << "Error={ ";
        for (const auto& e : solver.errors())
          std::cout << e << " ";
        std::cout << "} after " << solver.statistics().iterations << " iterations" << std::endl;
        auto evals = residual(x, g);
        for (size_t root = 0; root < solver.m_roots; root++) {
          std::cout << "Eigenvalue " << std::fixed << std::setprecision(9) << evals[root] << std::endl;
          //        solver.solution(root, x.front(), g.front(), Pcoeff.front());
          //        std::cout << "Eigenvector: (norm=" << std::sqrt(x[0].dot(x[0])) << "): ";
          //        for (size_t k = 0; k < n; k++)
          //          std::cout << " " << (x[0])[k];
          //        std::cout << std::endl;
        }
        EXPECT_THAT(solver.errors(),
                    ::testing::Pointwise(::testing::DoubleNear(solver.m_thresh), std::vector<double>(1, double(0))));
        EXPECT_THAT(evals, ::testing::Pointwise(::testing::DoubleNear(2e-9),
                                                std::vector<double>(expected_eigenvalues.data(),
                                                                    expected_eigenvalues.data() + solver.m_roots)));
      }
    }
  }
}
