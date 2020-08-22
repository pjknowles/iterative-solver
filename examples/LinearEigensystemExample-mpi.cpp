#include <molpro/linalg/IterativeSolver.h>
#include <fstream>
#include <iomanip>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <mpi.h>
#include <vector>
#include <cstring>

// Find lowest eigensolutions of a matrix obtained from an external file
using Rvector = molpro::linalg::array::DistrArrayMPI3;
//using Qvector = molpro::linalg::array::DistrArrayHDF5;
using Qvector = molpro::linalg::array::DistrArrayMPI3;
using Pvector = std::map<size_t, double>;
int n; // dimension of problem
int mpi_rank;
int mpi_size;
std::vector<double> hmat;

void action(size_t nwork, const std::vector<Rvector>& psc, std::vector<Rvector>& psg) {
  for (size_t k = 0; k < nwork; k++) {
    auto grange = psg[k].distribution().range(mpi_rank);
    auto gn = grange.second - grange.first;
    auto g_chunk = psg[k].local_buffer();
    for (auto gg=0; gg<gn; gg++) (*g_chunk)[gg]=0;
    for (int crank = 0; crank < mpi_size; crank++) {
      auto crange = psc[k].distribution().range(crank);
      auto cn = crange.second - crange.first;
      std::vector<double> c(cn);
      if (crank == mpi_rank)
        std::memcpy( c.data(), &(*psc[k].local_buffer())[0],cn* sizeof(double));
      MPI_Bcast(c.data(), cn, MPI_DOUBLE, crank, MPI_COMM_WORLD);
      for (size_t i = grange.first; i < grange.second; i++) {
        for (size_t j = crange.first; j < crange.second; j++)
          (*g_chunk)[i - grange.first] += hmat[j+i*n] * c[j - crange.first];
      }
    }
  }
}

void update(std::vector<Rvector>& psc, const std::vector<Rvector>& psg, size_t nwork,
            std::vector<double> shift = std::vector<double>()) {
  for (size_t k = 0; k < nwork; k++) {
    auto range = psg[k].distribution().range(mpi_rank);
    assert(range == psc[k].distribution().range(mpi_rank));
    auto c_chunk = psc[k].local_buffer();
    auto g_chunk = psg[k].local_buffer();
    for (size_t i = range.first; i < range.second; i++) {
      (*c_chunk)[i - range.first] -= (*g_chunk)[i - range.first] / (1e-12 - shift[k] + hmat[i+i*n]);
    }
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  for (const auto& file : std::vector<std::string>{"hf", "bh"}) {
    for (const auto& nroot : std::vector<int>{1, 2, 4}) {
      if (mpi_rank == 0) {
        std::string prefix{argv[0]};
        std::cout << prefix << std::endl;
        if (prefix.find_last_of("/") != std::string::npos)
          prefix.resize(prefix.find_last_of("/"));
        else
          prefix = ".";
        std::cout << prefix << std::endl;
        std::ifstream f(prefix + "/examples/" + file + ".hamiltonian");
        f >> n;
        molpro::cout << "\n*** " << file << " (dimension " << n << "), " << nroot << " roots, mpi_size = " << mpi_size
                     << std::endl;
        hmat.resize(n * n);
        for (auto i = 0; i < n * n; i++)
          f >> hmat[i];
      }
      MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (mpi_rank != 0)
        hmat.resize(n * n);
      MPI_Bcast(hmat.data(), n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      std::vector<double> diagonals;
      diagonals.reserve(n);
      for (auto i = 0; i < n; i++)
        diagonals.push_back(hmat[i+i*n]);
      molpro::linalg::LinearEigensystem<Rvector, Qvector, Pvector> solver;
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
      for (auto iter = 0; iter < 100; iter++) {
        action(nwork, x, g);
        nwork = solver.addVector(x, g, Pcoeff);
        if (mpi_rank == 0)
          solver.report();
        if (nwork == 0)
          break;
        update(x, g, nwork, solver.working_set_eigenvalues());
      }
      if (mpi_rank == 0) {
        std::cout << "Error={ ";
        for (const auto& e : solver.errors())
          std::cout << e << " ";
        std::cout << "} after " << solver.iterations() << " iterations" << std::endl;
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
          auto result = std::sqrt(handlers.rr().dot(g[root], g[root]));
          if (mpi_rank == 0)
            std::cout << " " << result;
        }
        if (mpi_rank == 0)
          std::cout << std::endl;
        if (mpi_rank == 0)
          std::cout << "Eigenvector orthonormality:\n";
        for (size_t root = 0; root < solver.m_roots; root++) {
          for (size_t soot = 0; soot < solver.m_roots; soot++) {
            auto result = handlers.rr().dot(x[root], x[soot]);
            if (mpi_rank == 0)
              std::cout << " " << result;
          }
          if (mpi_rank == 0)
            std::cout << std::endl;
          //        std::cout << "Eigenvector: (norm=" << std::sqrt(x[0].dot(x[0])) << "): ";
          //        for (size_t k = 0; k < n; k++)
          //          std::cout << " " << (x[0])[k];
          //        std::cout << std::endl;
        }
      }
    }
  }
  MPI_Finalize();
}
