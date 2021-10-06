#include "container.h"
#include <cstring>
#include <fstream>
#include <iomanip>
#include <molpro/linalg/itsolv/LinearEigensystemDavidson.h>
#include <vector>

// Find lowest eigensolutions of a matrix obtained from an external file
using Rvector = container<double>;

class ExampleProblem : public molpro::linalg::itsolv::Problem<Rvector> {
protected:
  std::vector<double> matrix;
  size_t n;

public:
  size_t size() const { return n; }
  using Problem::container_t;
  using Problem::value_t;
  ExampleProblem(std::string filename) {
    std::ifstream f(filename);
    f >> n;
    //    molpro::cout << "\n*** " << filename << " (dimension " << n << "), " << nroot << " roots" << std::endl;
    matrix.resize(n * n);
    for (size_t i = 0; i < n * n; i++)
      f >> matrix[i];
  }

  bool diagonals(container_t &d) const override {
    for (size_t i = 0; i < n; i++)
      d.data()[i] = matrix[i * (n + 1)];
    return true;
  }

  void precondition(const VecRef<container_t> &action, const std::vector<double> &shift,
                    const container_t &diagonals) const override {
    for (size_t k = 0; k < action.size(); k++) {
      auto &a = action[k].get();
      for (size_t i = 0; i < n; i++)
        a.data()[i] /= (diagonals.data()[i] - shift[k] + 1e-15);
    }
  }

  void action(const CVecRef<container_t> &parameters, const VecRef<container_t> &actions) const override {
    for (size_t k = 0; k < parameters.size(); k++) {
      const auto &v = parameters[k].get().data();
      auto &a = actions[k].get().data();
      for (size_t i = 0; i < n; i++) {
        a[i] = 0;
        for (size_t j = 0; j < n; j++)
          a[i] += matrix[j + n * i] * v[j];
      }
    }
  }
};

int main(int argc, char *argv[]) {
  for (const auto &file : std::vector<std::string>{"hf", "bh"}) {
    for (const auto &nroot : std::vector<int>{1, 2}) {
      std::string prefix{argv[0]};
      if (prefix.find_last_of("/") != std::string::npos)
        prefix.resize(prefix.find_last_of("/"));
      else
        prefix = ".";
      auto fullfile = prefix + "/../" + file + ".hamiltonian";
      ExampleProblem problem(fullfile);
      molpro::cout << "\n*** " << fullfile << " (dimension " << problem.size() << "), " << nroot << " roots"
                   << std::endl;
      molpro::linalg::itsolv::LinearEigensystemDavidson<Rvector> itsolver;
      itsolver.set_n_roots(nroot);
      itsolver.set_convergence_threshold(1.0e-12);
      std::vector<Rvector> g;
      std::vector<Rvector> x;
      //      auto low_diagonals =
      for (size_t root = 0; root < itsolver.n_roots(); root++) {
        x.emplace_back(problem.size());
        g.emplace_back(problem.size());
      }
      itsolver.solve(x, g, problem, true);
      {
        std::cout << "Error={ ";
        for (const auto &e : itsolver.errors())
          std::cout << e << " ";
        std::cout << "} after " << itsolver.statistics().iterations << " iterations" << std::endl;
        std::cout << "Statistics: " << itsolver.statistics() << std::endl;
        for (size_t root = 0; root < itsolver.n_roots(); root++) {
          std::cout << "Eigenvalue " << std::fixed << std::setprecision(9) << itsolver.eigenvalues()[root] << std::endl;
        }
      }
      {
        auto working_set = itsolver.working_set();
        working_set.resize(itsolver.n_roots());
        std::iota(working_set.begin(), working_set.end(), 0);
        itsolver.solution(working_set, x, g);
        std::cout << "Residual norms:";
        for (size_t root = 0; root < itsolver.n_roots(); root++)
          std::cout << " " << std::sqrt(g[root].dot(g[root]));
        std::cout << std::endl;
        std::cout << "Eigenvector orthonormality:\n";
        for (size_t root = 0; root < itsolver.n_roots(); root++) {
          for (size_t soot = 0; soot < itsolver.n_roots(); soot++)
            std::cout << " " << x[root].dot(x[soot]);
          std::cout << std::endl;
        }
      }
    }
  }
}