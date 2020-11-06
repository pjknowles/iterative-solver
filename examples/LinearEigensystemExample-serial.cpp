#include <cstring>
#include <fstream>
#include <iomanip>
#include <molpro/linalg/IterativeSolver.h>
#include <molpro/linalg/itsolv/LinearEigensystemA.h>
#include <vector>

// Find lowest eigensolutions of a matrix obtained from an external file
using Rvector = std::vector<double>;
using Qvector = std::vector<double>;
using Pvector = std::map<size_t, double>;
int n; // dimension of problem
std::vector<double> hmat;

void action(size_t nwork, const std::vector<Rvector>& psc, std::vector<Rvector>& psg) {
  for (size_t k = 0; k < nwork; k++) {
    std::fill(psg[k].begin(), psg[k].end(), 0);
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
        psg[k][j] += hmat[j + i * n] * psc[k][i];
  }
}

void update(std::vector<Rvector>& psc, std::vector<Rvector>& psg, size_t nwork,
            std::vector<double> shift = std::vector<double>()) {
  for (size_t k = 0; k < nwork; k++) {
    for (size_t i = 0; i < n; i++)
      psg[k][i] *= 1 / (1e-12 - shift[k] + hmat[i + i * n]);
  }
}

int main(int argc, char* argv[]) {
  for (const auto& file : std::vector<std::string>{"bh"}) {
    for (const auto& nroot : std::vector<int>{4}) {
      std::string prefix{argv[0]};
      if (prefix.find_last_of("/") != std::string::npos)
        prefix.resize(prefix.find_last_of("/"));
      else
        prefix = ".";
      std::ifstream f(prefix + "/examples/" + file + ".hamiltonian");
      f >> n;
      if (nroot > n)
        continue;
      molpro::cout << "\n*** " << file << " (dimension " << n << "), " << nroot << " roots" << std::endl;
      hmat.resize(n * n);
      for (auto i = 0; i < n * n; i++)
        f >> hmat[i];
      std::vector<double> diagonals;
      diagonals.reserve(n);
      for (auto i = 0; i < n; i++)
        diagonals.push_back(hmat[i + i * n]);
      auto handlers = std::make_shared<molpro::linalg::itsolv::ArrayHandlers<Rvector, Qvector, Pvector>>();
      molpro::linalg::itsolv::LinearEigensystemA<Rvector, Qvector, Pvector> solver(handlers);
      solver.set_n_roots(nroot);
      solver.set_convergence_threshold(1.0e-12);
      solver.propose_rspace_norm_thresh = 1.0e-14;
      solver.max_size_qspace = 10;
      solver.set_reset_D(50);
      solver.logger->max_trace_level = molpro::linalg::itsolv::Logger::Info;
      solver.logger->max_warn_level = molpro::linalg::itsolv::Logger::Error;
      solver.logger->data_dump = true;
      std::vector<Rvector> g;
      std::vector<Rvector> x;
      for (size_t root = 0; root < solver.n_roots(); root++) {
        x.emplace_back(n);
        g.emplace_back(n);
        std::fill(x.back().begin(), x.back().end(), 0);
        auto guess = std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin(); // initial guess
        x.back()[guess] = 1;
        *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
      }
      std::vector<std::vector<double>> Pcoeff(solver.n_roots());
      int nwork = solver.n_roots();
      bool done = false;
      for (auto iter = 0; iter < 10 && !done; iter++) {
        if (false) {
          solver.logger->max_trace_level = molpro::linalg::itsolv::Logger::Info;
          solver.logger->max_warn_level = molpro::linalg::itsolv::Logger::Error;
          solver.logger->data_dump = true;
        }
        action(nwork, x, g);
        nwork = solver.add_vector(x, g);
        solver.report();
        done = nwork == 0;
        if (nwork != 0) {
          update(x, g, nwork, solver.working_set_eigenvalues());
          nwork = solver.end_iteration(x, g);
        }
      }
      {
        std::cout << "Error={ ";
        for (const auto& e : solver.errors())
          std::cout << e << " ";
        std::cout << "} after " << solver.statistics().iterations << " iterations" << std::endl;
        std::cout << "Statistics: " << solver.statistics() << std::endl;
        for (size_t root = 0; root < solver.n_roots(); root++) {
          std::cout << "Eigenvalue " << std::fixed << std::setprecision(9) << solver.eigenvalues()[root] << std::endl;
        }
      }
      {
        auto working_set = solver.working_set();
        working_set.resize(solver.n_roots());
        std::iota(working_set.begin(), working_set.end(), 0);
        solver.solution(working_set, x, g);
        std::cout << "Residual norms:";
        for (size_t root = 0; root < solver.n_roots(); root++)
          std::cout << " " << std::sqrt(handlers->rr().dot(g[root], g[root]));
        std::cout << std::endl;
        std::cout << "Eigenvector orthonormality:\n";
        for (size_t root = 0; root < solver.n_roots(); root++) {
          for (size_t soot = 0; soot < solver.n_roots(); soot++)
            std::cout << " " << handlers->rr().dot(x[root], x[soot]);
          std::cout << std::endl;
        }
      }
    }
  }
}
