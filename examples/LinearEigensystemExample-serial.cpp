#include <cstring>
#include <fstream>
#include <iomanip>
#include <molpro/linalg/IterativeSolver.h>
#include <vector>

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  bool init{true};
  for (const auto& element : v) {
    o << (init ? "{" : " ") << element;
    init = false;
  }
  o << "}";
  return o;
}

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
        psg[k][i] += hmat[i + j * n] * psc[k][j];
  }
}

void actionP(size_t nwork, const std::vector<Pvector>& pspace, const std::vector<std::vector<double>> Pcoeff, std::vector<Rvector>& psg) {
  for (size_t k = 0; k < nwork; k++) {
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < pspace.size(); j++)
        psg[k][i] += hmat[i + pspace[j].begin()->first * n] * Pcoeff[k][j];
  }
}

void update(std::vector<Rvector>& psc, const std::vector<Rvector>& psg, size_t nwork,
            std::vector<double> shift = std::vector<double>()) {
  for (size_t k = 0; k < nwork; k++) {
    for (size_t i = 0; i < n; i++)
      psc[k][i] -= psg[k][i] / (1e-12 - shift[k] + hmat[i + i * n]);
  }
}

int main(int argc, char* argv[]) {
  for (const auto& file : std::vector<std::string>{"hf", "bh"}) {
    for (const auto& nroot : std::vector<int>{1, 2, 4}) {
      std::string prefix{argv[0]};
      std::cout << prefix << std::endl;
      if (prefix.find_last_of("/") != std::string::npos)
        prefix.resize(prefix.find_last_of("/"));
      else
        prefix = ".";
      std::cout << prefix << std::endl;
      std::ifstream f(prefix + "/examples/" + file + ".hamiltonian");
      f >> n;
      molpro::cout << "\n*** " << file << " (dimension " << n << "), " << nroot << " roots" << std::endl;
      hmat.resize(n * n);
      for (auto i = 0; i < n * n; i++)
        f >> hmat[i];
      //      std::cout << "hmat "<<hmat<<std::endl;
      //      for (const auto& nP : std::vector<size_t>{0, 1, 3, 5}) {
      for (const auto& nP : std::vector<size_t>{1}) {
        std::vector<double> diagonals;
        diagonals.reserve(n);
        for (auto i = 0; i < n; i++)
          diagonals.push_back(hmat[i + i * n]);
        //      std::cout << "diagonal matrix elements "<<diagonals<<std::endl;
        if (nP > n)
          continue;
        std::cout << "P-space dimension = " << nP << std::endl;
        auto handlers = std::make_shared<molpro::linalg::iterativesolver::ArrayHandlers<Rvector, Qvector, Pvector>>();
        auto solver = molpro::linalg::LinearEigensystem<Rvector, Qvector, Pvector>(handlers);
        solver.m_verbosity = 3;
        solver.m_roots = nroot;
        solver.m_thresh = 1e-9;
        std::vector<Rvector> g;
        std::vector<Rvector> x;
        std::vector<std::vector<double>> Pcoeff;
        double d0 = *std::min_element(diagonals.begin(), diagonals.end());
        for (size_t root = 0; root < solver.m_roots; root++) {
          x.emplace_back(n);
          g.emplace_back(n);
          std::fill(x.back().begin(), x.back().end(), 0);
          auto guess = std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin(); // initial guess
          x.back()[guess] = 1;
          *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
          Pcoeff.emplace_back(nP);
        }
        for (auto i = 0; i < n; i++)
          diagonals[i] = 1 / (1e-12 + hmat[i + i * n] - d0);
        std::cout << "diagonal matrix elements " << diagonals << std::endl;
        auto selection =
            handlers->rr().select_max_dot( nP, diagonals, diagonals);
//        molpro::linalg::array::util::select_max_dot<std::vector<double>, std::vector<double>, double, double>( nP, diagonals, diagonals);
        std::vector<Pvector> pspace;
        std::vector<double> hpp;
        hpp.reserve(nP * nP);
        size_t p = 0;
        for (const auto& select : selection) {
          pspace.emplace_back();
          pspace.back()[select.first] = 1;
          std::cout << "P space element " << select.first << std::endl;
          for (const auto& select2 : selection) {
            hpp.push_back(hmat[select2.first + n * select.first]);
          }
        }
        std::cout << "hpp" << hpp << std::endl;
        int nwork = solver.m_roots;
        for (auto iter = 0; iter < 10; iter++) {
          if (iter == 0 && nP > 0) {
            nwork = solver.addP(pspace, hpp.data(), x, g, Pcoeff);
          } else {
            action(nwork, x, g);
            std::cout << "before addVector x=" << x << std::endl;
            std::cout << "before addVector g=" << g << std::endl;
            nwork = solver.addVector(x, g, Pcoeff);
          }
          std::cout << "nwork=" << nwork << std::endl;
          std::cout << "after add* x=" << x << std::endl;
          std::cout << "after add* g=" << g << std::endl;
          for (const auto& Pcoefff : Pcoeff)
            std::cout << "after add* Pcoeff column " << Pcoefff << std::endl;
          solver.report();
          if (nwork == 0)
            break;
          actionP(nwork,pspace,Pcoeff,g);
          update(x, g, nwork, solver.working_set_eigenvalues());
        }
        {
          std::cout << "Error={ ";
          for (const auto& e : solver.errors())
            std::cout << e << " ";
          std::cout << "} after " << solver.statistics().iterations << " iterations" << std::endl;
          for (size_t root = 0; root < solver.m_roots; root++) {
            std::cout << "Eigenvalue " << std::fixed << std::setprecision(9) << solver.eigenvalues()[root] << std::endl;
          }
        }
        {
          auto working_set = solver.working_set();
          working_set.resize(solver.m_roots);
          std::iota(working_set.begin(), working_set.end(), 0);
          solver.solution(working_set, x, g, Pcoeff);
          std::cout << "Residual norms:";
          for (size_t root = 0; root < solver.m_roots; root++)
            std::cout << " " << std::sqrt(handlers->rr().dot(g[root], g[root]));
          std::cout << std::endl;
          std::cout << "Eigenvector orthonormality:\n";
          for (size_t root = 0; root < solver.m_roots; root++) {
            for (size_t soot = 0; soot < solver.m_roots; soot++)
              std::cout << " " << handlers->rr().dot(x[root], x[soot]);
            std::cout << std::endl;
          }
        }
      }
    }
  }
}
