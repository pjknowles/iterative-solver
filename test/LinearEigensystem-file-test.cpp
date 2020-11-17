#include "test.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "molpro/linalg/IterativeSolver.h"
#include "molpro/linalg/PagedArray.h"
#include "molpro/linalg/SimpleArray.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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
#include <molpro/linalg/itsolv/LinearEigensystemA.h>

using molpro::linalg::array::ArrayHandler;
using molpro::linalg::array::ArrayHandlerIterable;
using molpro::linalg::array::ArrayHandlerIterableSparse;
using molpro::linalg::array::ArrayHandlerSparse;
using molpro::linalg::itsolv::ArrayHandlers;
using Rvector = std::vector<double>;
using Qvector = std::vector<double>;
using Pvector = std::map<size_t, double>;
// Find lowest eigensolution of a matrix obtained from an external file
// Storage of vectors in-memory via class Rvector
using scalar = double;
size_t n;
std::vector<double> hmat;
std::vector<double> expected_eigenvalues;

void load_matrix(int dimension, const std::string& type = "", double param = 1) {
  n = dimension;
  hmat.resize(n * n);
  hmat.assign(n * n, 1);
  for (int i = 0; i < n; i++)
    hmat[i * (n + 1)] = i * param;
}
void load_matrix(const std::string& file) {
  std::ifstream f(std::string{"./"} + file + ".hamiltonian");
  f >> n;
  hmat.resize(n * n);
  for (auto i = 0; i < n * n; i++)
    f >> hmat[i];
}

using pv = Rvector;

scalar matrix(const size_t i, const size_t j) { return hmat[i * n + j]; }

void action(const std::vector<Rvector>& psx, std::vector<Rvector>& outputs) {
  for (size_t k = 0; k < psx.size(); k++) {
    for (size_t i = 0; i < n; i++) {
      outputs[k][i] = 0;
      for (size_t j = 0; j < n; j++)
        outputs[k][i] += matrix(i, j) * psx[k][j];
    }
  }
}
template <class scalar>
scalar dot(const std::vector<scalar>& a, const std::vector<scalar>& b) {
  return std::inner_product(std::begin(a), std::end(a), std::begin(b), scalar(0));
}
std::vector<double> residual(const std::vector<Rvector>& psx, std::vector<Rvector>& outputs) {

  std::vector<double> evals;
  for (size_t k = 0; k < psx.size(); k++) {
    double eval = 0;
    for (size_t i = 0; i < n; i++) {
      outputs[k][i] = 0;
      for (size_t j = 0; j < n; j++)
        outputs[k][i] += matrix(i, j) * psx[k][j];
      eval += outputs[k][i] * psx[k][i];
    }
    eval /= dot(psx[k], psx[k]);
    for (size_t i = 0; i < n; i++) {
      outputs[k][i] -= eval * psx[k][i];
    }
    evals.push_back(eval);
  }
  return evals;
}

void update(std::vector<Rvector>& psg, std::vector<scalar> shift = std::vector<scalar>()) {
  //  for (size_t k = 0; k < psc.size(); k++) {
  //   molpro::cout << "update root "<<k<<", shift="<<shift[k]<<": ";
  //    for (size_t i = 0; i < n; i++)
  //      molpro::cout <<" "<< -psg[k][i] / (1e-12-shift[k] + matrix(i,i));
  //    molpro::cout<<std::endl;
  //  }
  for (size_t k = 0; k < psg.size() && k < shift.size(); k++)
    for (size_t i = 0; i < n; i++)
      psg[k][i] = -psg[k][i] / (1e-12 - shift[k] + matrix(i, i));
}

void test_eigen(const std::string& title = "") {
  {
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> he(hmat.data(), n, n);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> esolver(he);
    auto ev = esolver.eigenvalues();
    //              std::cout << "actual eigenvalues" << ev << std::endl;
    expected_eigenvalues.clear();
    for (int i = 0; i < n; i++)
      expected_eigenvalues.push_back(ev[i]);
  }
  for (int nroot = 1; nroot <= std::min(n, n / 2 + 3) && nroot <= 23; nroot++) {
    for (auto np = 0; np <= 0; np += 4) {
      molpro::cout << "\n\n*** " << title << ", " << nroot << " roots, problem dimension " << n << ", pspace dimension "
                   << np << std::endl;

      auto handlers = std::make_shared<molpro::linalg::itsolv::ArrayHandlers<Rvector, Qvector, Pvector>>();
      molpro::linalg::itsolv::LinearEigensystemA<Rvector, Qvector, Pvector> solver(handlers);
      solver.set_n_roots(nroot);
      solver.set_convergence_threshold(1.0e-10);
      solver.propose_rspace_norm_thresh = 1.0e-14;
      solver.propose_rspace_svd_thresh = 1.0e-4;
      solver.max_size_qspace = std::min(int(n), std::min(1000, 3 * nroot));
      solver.set_reset_D(10);
      molpro::cout << "convergence threshold = " << solver.convergence_threshold() << ", svd thresh"
                   << solver.propose_rspace_svd_thresh << ", norm thresh" << solver.propose_rspace_norm_thresh
                   << ", max size of Q = " << solver.max_size_qspace << ", reset D = " << solver.get_reset_D()
                   << std::endl;
      solver.logger->max_trace_level = molpro::linalg::itsolv::Logger::None;
      solver.logger->max_warn_level = molpro::linalg::itsolv::Logger::Error;
      solver.logger->data_dump = false;
      std::vector<Rvector> g;
      std::vector<Rvector> x;

      std::vector<size_t> guess;
      {
        std::vector<double> diagonals;
        for (auto i = 0; i < n; i++)
          diagonals.push_back(matrix(i, i));
        for (size_t root = 0; root < solver.n_roots(); root++) {
          x.emplace_back(n, 0);
          g.emplace_back(n);
          guess.push_back(std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin()); // initial guess
          *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
          x.back()[guess.back()] = 1; // initial guess
        }
      }
      //        std::cout << "guess: " << guess << std::endl;
      int nwork = solver.n_roots();

      std::vector<std::vector<scalar>> Pcoeff(solver.n_roots());
      std::vector<scalar> PP;
      std::vector<Pvector> pspace;

      for (auto iter = 0; iter < 100; iter++) {
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
          using vectorP = std::vector<scalar>;
          using molpro::linalg::itsolv::CVecRef;
          using molpro::linalg::itsolv::VecRef;
//          using fapply_on_p_type = std::function<void(const std::vector<VectorP>&, const CVecRef<P>&, const VecRef<R>&)>;
          std::function<void(const std::vector<vectorP>&, const CVecRef<Pvector>&, const VecRef<Rvector>&)>
//          auto
          apply_p = [](const std::vector<vectorP>& pvectors, const CVecRef<Pvector>& pspace, VecRef<Rvector>& action) {
            for (size_t i = 0; i < pvectors.size(); i++) {
              auto& actioni = (action[i]).get();
              for (size_t pi = 0; pi < pspace.size(); pi++) {
                const auto& p = pspace[pi].get();
                for (const auto& pel : p)
                  for (size_t j = 0; j < n; j++)
                    actioni[j] += hmat[j + n*pel.first] * pel.second * pvectors[i][pi];
              }
            }
          };
          nwork = solver.add_p(molpro::linalg::itsolv::cwrap(pspace),
                               molpro::linalg::array::span::Span<double>(PP.data(), PP.size()),
                               molpro::linalg::itsolv::wrap(x), molpro::linalg::itsolv::wrap(g), Pcoeff, apply_p);
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
          nwork = solver.add_vector(x, g, Pcoeff);
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
          //            solver.m_errors[solver.working_set()[k]] = std::sqrt(g[k].dot(g[k]));
          //            molpro::cout << "; new error " << solver.m_errors[solver.working_set()[k]] << std::endl;
        }
        //          for (auto root = 0; root < nwork; root++) {
        //            std::cout << "after p gradient g:";
        //            for (auto i = 0; i < n; i++)
        //              std::cout << " " << g[root][i];
        //            std::cout << std::endl;
        //          }
        update(g, solver.working_set_eigenvalues());
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
        //          if (*std::max_element(solver.errors().begin(), solver.errors().end()) < solver.m_thresh)
        if (solver.end_iteration(x, g) == 0)
          break;
      }
      std::cout << "Error={ ";
      for (const auto& e : solver.errors())
        std::cout << e << " ";
      std::cout << "} after " << solver.statistics().iterations << " iterations, " << solver.statistics().r_creations
                << " R vectors" << std::endl;
      //        for (size_t root = 0; root < solver.n_roots(); root++) {
      //          std::cout << "Eigenvalue " << std::fixed << std::setprecision(9) << solver.eigenvalues()[root] <<
      //          std::endl;
      //        solver.solution(root, x.front(), g.front(), Pcoeff.front());
      //        std::cout << "Eigenvector: (norm=" << std::sqrt(x[0].dot(x[0])) << "): ";
      //        for (size_t k = 0; k < n; k++)
      //          std::cout << " " << (x[0])[k];
      //        std::cout << std::endl;
      //        }
      EXPECT_THAT(solver.errors(), ::testing::Pointwise(::testing::DoubleNear(2 * solver.convergence_threshold()),
                                                        std::vector<double>(nroot, double(0))));
      std::cout << "expected eigenvalues "
                << std::vector<double>(expected_eigenvalues.data(), expected_eigenvalues.data() + solver.n_roots())
                << std::endl;
      std::cout << "obtained eigenvalues " << solver.eigenvalues() << std::endl;
      EXPECT_THAT(solver.eigenvalues(),
                  ::testing::Pointwise(::testing::DoubleNear(2e-9),
                                       std::vector<double>(expected_eigenvalues.data(),
                                                           expected_eigenvalues.data() + solver.n_roots())));
      EXPECT_LE(solver.statistics().r_creations, (nroot + 1) * 15);
    }
  }
}

TEST(IterativeSolver, file_eigen) {
  for (const auto& file : std::vector<std::string>{"bh"}) {
    load_matrix(file);
    test_eigen(file);
  }
}

TEST(IterativeSolver, n_eigen) {
  size_t n = 1000;
  double param = 1;
  //  for (auto param : std::vector<double>{.01, .1, 1, 10, 100}) {
  //  for (auto param : std::vector<double>{.01, .1, 1}) {
  for (auto param : std::vector<double>{1}) {
    load_matrix(n, "", param);
    test_eigen(std::to_string(n) + "/" + std::to_string(param));
  }
}

TEST(IterativeSolver, small_eigen) {
  for (int n = 1; n < 20; n++) {
    double param = 1;
    load_matrix(n, "", param);
    test_eigen(std::to_string(n) + "/" + std::to_string(param));
  }
}

TEST(IterativeSolver, symmetry_eigen) {
  for (int n = 1; n < 6; n++) {
    double param = 1;
    load_matrix(n, "", param);
    for (auto i = 0; i < n; i++)
      for (auto j = 0; j < n; j++)
        if (((i % 3) == 0 and (j % 3) != 0) or ((j % 3) == 0 and (i % 3) != 0))
          hmat[j + i * n] = 0;
    if (false) {
      std::cout << "hmat " << std::endl;
      for (auto i = 0; i < n; i++) {
        std::cout << "i%3" << i % 3 << std::endl;
        for (auto j = 0; j < n; j++)
          std::cout << " " << hmat[j + i * n];
        std::cout << std::endl;
      }
    }
    test_eigen(std::to_string(n) + "/sym/" + std::to_string(param));
  }
}

TEST(IterativeSolver, DISABLED_file_optimize_eigenvalue) {
  for (const auto& file : std::vector<std::string>{"bh"}) {
    load_matrix(file);
    {
      molpro::cout << "\n\n*** " << file << ", BFGS for lowest eigensolution, problem dimension " << n << std::endl;
      auto rr = std::make_shared<ArrayHandlerIterable<Rvector>>();
      auto qq = std::make_shared<ArrayHandlerIterable<Rvector>>();
      auto pp = std::make_shared<ArrayHandlerSparse<std::map<size_t, double>>>();
      auto rq = std::make_shared<ArrayHandlerIterable<Rvector>>();
      auto rp = std::make_shared<ArrayHandlerIterableSparse<Rvector, std::map<size_t, double>>>();
      auto qr = std::make_shared<ArrayHandlerIterable<Rvector>>();
      auto qp = std::make_shared<ArrayHandlerIterableSparse<Rvector, std::map<size_t, double>>>();
      auto handlers =
          std::make_shared<ArrayHandlers<Rvector, Rvector, std::map<size_t, double>>>(rr, qq, pp, rq, rp, qr, qp);
      auto solver = molpro::linalg::Optimize<Rvector>{handlers};
      solver.m_verbosity = 1;
      solver.m_thresh = 1e-9;
      std::vector<Rvector> g;
      std::vector<Rvector> x;
      std::vector<size_t> guess;
      {
        std::vector<double> diagonals;
        for (auto i = 0; i < n; i++)
          diagonals.push_back(matrix(i, i));
        for (size_t root = 0; root < solver.m_roots; root++) {
          x.emplace_back(n);
          g.emplace_back(n);
          x.back().assign(x.back().size(), 0);
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
        update(g, {eval});
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

TEST(IterativeSolver, DISABLED_file_diis_eigenvalue) {
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
        auto rr = std::make_shared<ArrayHandlerIterable<Rvector>>();
        auto qq = std::make_shared<ArrayHandlerIterable<Rvector>>();
        auto pp = std::make_shared<ArrayHandlerSparse<std::map<size_t, double>>>();
        auto rq = std::make_shared<ArrayHandlerIterable<Rvector>>();
        auto rp = std::make_shared<ArrayHandlerIterableSparse<Rvector, std::map<size_t, double>>>();
        auto qr = std::make_shared<ArrayHandlerIterable<Rvector>>();
        auto qp = std::make_shared<ArrayHandlerIterableSparse<Rvector, std::map<size_t, double>>>();
        auto handlers =
            std::make_shared<ArrayHandlers<Rvector, Rvector, std::map<size_t, double>>>(rr, qq, pp, rq, rp, qr, qp);
        auto solver = molpro::linalg::DIIS<Rvector>{handlers};
        solver.m_verbosity = 1;
        solver.m_thresh = 1e-9;
        std::vector<Rvector> g;
        std::vector<Rvector> x;
        std::vector<size_t> guess;
        {
          std::vector<double> diagonals;
          for (auto i = 0; i < n; i++)
            diagonals.push_back(matrix(i, i));
          for (size_t root = 0; root < solver.m_roots; root++) {
            x.emplace_back(n);
            g.emplace_back(n);
            x.back().assign(x.back().size(), 0);
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
          update(g, {eval});
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
