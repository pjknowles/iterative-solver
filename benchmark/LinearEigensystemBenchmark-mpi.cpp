#include <cstring>
#include <fstream>
#include <iomanip>
#include <molpro/linalg/IterativeSolver.h>
#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <mpi.h>
#include <vector>

// Find lowest eigensolutions of M(i,j) = alpha*(i+1)*delta(i,j) + 1
using scalar = double;
size_t n;                   // dimension of problem
constexpr scalar alpha = 1; // separation of diagonal elements
constexpr bool collective_algorithm = false;
// scalar matrix(const size_t i, const size_t j) { return (i == j ? alpha * (i + 1) : 0) + (i + j)%((n+7)/6); }
scalar matrix(const size_t i, const size_t j) { return (i == j ? alpha * (i + 1) : 0) + 1; }
using Rvector = molpro::linalg::array::DistrArrayMPI3;
// using Qvector = molpro::linalg::array::DistrArrayHDF5;
using Qvector = molpro::linalg::array::DistrArrayMPI3;
using Pvector = std::map<size_t, double>;
int mpi_rank;
int mpi_size;
std::vector<double> hmat;

void action(size_t nwork, const std::vector<Rvector>& psc, std::vector<Rvector>& psg) {
  if (not collective_algorithm)
    MPI_Barrier(MPI_COMM_WORLD);
  for (size_t k = 0; k < nwork; k++) {
    auto grange = psg[k].distribution().range(mpi_rank);
    auto gn = grange.second - grange.first;
    auto g_chunk = psg[k].local_buffer();
    for (auto gg = 0; gg < gn; gg++)
      (*g_chunk)[gg] = 0;
    for (int cseg = 0; cseg < mpi_size; cseg++) {
      auto crank = collective_algorithm ? cseg : (mpi_size + cseg - mpi_rank) % mpi_size;
      auto crange = psc[k].distribution().range(crank);
      auto cn = crange.second - crange.first;
      std::vector<double> c(cn);
      if (crank == mpi_rank or (not collective_algorithm))
        psc[k].get(crange.first, crange.second, c.data());
      if (collective_algorithm)
        MPI_Bcast(c.data(), cn, MPI_DOUBLE, crank, MPI_COMM_WORLD);
      for (size_t i = grange.first; i < grange.second; i++) {
        for (size_t j = crange.first; j < crange.second; j++)
          (*g_chunk)[i - grange.first] += matrix(j, i) * c[j - crange.first];
      }
    }
  }
}

void update(std::vector<Rvector>& psc, const std::vector<Rvector>& psg, size_t nwork,
            std::vector<double> shift = std::vector<double>()) {
  for (size_t k = 0; k < nwork; k++) {
    auto range = psg[k].distribution().range(mpi_rank);
    assert(psg[k].compatible(psc[k]));
    auto c_chunk = psc[k].local_buffer();
    const auto g_chunk = psg[k].local_buffer();
    for (size_t i = range.first; i < range.second; i++) {
      (*c_chunk)[i - range.first] -= (*g_chunk)[i - range.first] / (1e-12 - shift[k] + matrix(i, i));
    }
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  for (const auto& nn : std::vector<int>{1, 2, 4, 10, 100, 1000, 10000}) {
    n = nn;
    for (const auto& nroot : std::vector<int>{1, 2, 4}) {
      if (nroot > n)
        break;
      if (mpi_rank == 0) {
        molpro::cout << "\n*** dimension " << n << ", " << nroot << " roots, mpi_size = " << mpi_size << std::endl;
      }
      std::vector<double> diagonals;
      diagonals.reserve(n);
      for (auto i = 0; i < n; i++)
        diagonals.push_back(matrix(i, i));
      auto solver = molpro::linalg::LinearEigensystem<Rvector, Qvector, Pvector>(
          std::make_shared<molpro::linalg::iterativesolver::ArrayHandlers<Rvector, Qvector, Pvector>>());
      auto handlers = solver.handlers();
      solver.m_verbosity = 1;
      solver.m_roots = nroot;
      solver.m_thresh = 1e-9;
      std::vector<Rvector> g;
      std::vector<Rvector> x;
      for (size_t root = 0; root < solver.m_roots; root++) {
        x.emplace_back(n, MPI_COMM_WORLD);
        g.emplace_back(n, MPI_COMM_WORLD);
        x.back().allocate_buffer();
        g.back().allocate_buffer();
        x.back().fill(0);
        auto guess = std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin(); // initial guess
        if (x.back().distribution().cover(guess) < n)
          x.back().set(guess, 1);
        *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
      }
      std::vector<std::vector<double>> Pcoeff(solver.m_roots);
      int nwork = solver.m_roots;
      for (auto iter = 0; iter < 100 && nwork > 0; iter++) {
        action(nwork, x, g);
        nwork = solver.addVector(x, g, Pcoeff);
        if (mpi_rank == 0)
          solver.report();
        update(x, g, nwork, solver.working_set_eigenvalues());
      }
      if (mpi_rank == 0) {
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
        solver.solution(working_set, x, g);
        if (mpi_rank == 0)
          std::cout << "Residual norms:";
        for (size_t root = 0; root < solver.m_roots; root++) {
          auto result = std::sqrt(handlers->rr().dot(g[root], g[root]));
          if (mpi_rank == 0)
            std::cout << " " << result;
        }
        if (mpi_rank == 0)
          std::cout << std::endl;
        if (mpi_rank == 0)
          std::cout << "Eigenvector orthonormality:\n";
        for (size_t root = 0; root < solver.m_roots; root++) {
          for (size_t soot = 0; soot < solver.m_roots; soot++) {
            auto result = handlers->rr().dot(x[root], x[soot]);
            if (mpi_rank == 0)
              std::cout << " " << result;
          }
          if (mpi_rank == 0)
            std::cout << std::endl;
        }
      }
      if (mpi_rank == 0)
        std::cout << solver.statistics() << std::endl;
    }
  }
  MPI_Finalize();
}
