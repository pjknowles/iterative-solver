#include <cstring>
#include <fstream>
#include <iomanip>
#include <molpro/linalg/IterativeSolver.h>
#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/itsolv/LinearEigensystemA.h>
#include <mpi.h>
#include <vector>

// Find lowest eigensolutions of a matrix obtained from an external file
using Rvector = molpro::linalg::array::DistrArrayMPI3;
using Qvector = molpro::linalg::array::DistrArrayHDF5;
using Pvector = std::map<size_t, double>;
int n; // dimension of problem
constexpr bool collective_algorithm = false;
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
          (*g_chunk)[i - grange.first] += hmat[j + i * n] * c[j - crange.first];
      }
    }
  }
}

void update(std::vector<Rvector>& psg, size_t nwork, std::vector<double> shift = std::vector<double>()) {
  for (size_t k = 0; k < nwork; k++) {
    auto range = psg[k].distribution().range(mpi_rank);
    auto g_chunk = psg[k].local_buffer();
    for (size_t i = range.first; i < range.second; i++) {
      (*g_chunk)[i - range.first] *= 1 / (1e-12 - shift[k] + hmat[i + i * n]);
    }
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  for (const auto& file : std::vector<std::string>{"hf", "bh"}) {
    for (const auto& nroot : std::vector<int>{1, 2, 4, 5}) {
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
        diagonals.push_back(hmat[i + i * n]);
      auto handlers = std::make_shared<molpro::linalg::itsolv::ArrayHandlers<Rvector, Qvector, Pvector>>();
      molpro::linalg::itsolv::LinearEigensystemA<Rvector, Qvector, Pvector> solver(handlers);
      solver.set_n_roots(nroot);
      solver.set_convergence_threshold(1.0e-12);
      solver.propose_rspace_norm_thresh = 1.0e-14;
      solver.set_max_size_qspace(10);
      solver.set_reset_D(50);
      solver.logger->max_trace_level = molpro::linalg::itsolv::Logger::None;
      solver.logger->max_warn_level = molpro::linalg::itsolv::Logger::Error;
      solver.logger->data_dump = false;
      std::vector<Rvector> g;
      std::vector<Rvector> x;
      for (size_t root = 0; root < solver.n_roots(); root++) {
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
      std::vector<std::vector<double>> Pcoeff(solver.n_roots());
      int nwork = solver.n_roots();
      bool done = false;
      for (auto iter = 0; iter < 100 && !done; iter++) {
        action(nwork, x, g);
        nwork = solver.add_vector(x, g);
        if (mpi_rank == 0)
          solver.report();
        done = nwork == 0;
        if (nwork != 0) {
          update(g, nwork, solver.working_set_eigenvalues());
          nwork = solver.end_iteration(x, g);
        }
      }
      if (mpi_rank == 0) {
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
        if (mpi_rank == 0)
          std::cout << "Residual norms:";
        for (size_t root = 0; root < solver.n_roots(); root++) {
          auto result = std::sqrt(handlers->rr().dot(g[root], g[root]));
          if (mpi_rank == 0)
            std::cout << " " << result;
        }
        if (mpi_rank == 0)
          std::cout << std::endl;
        if (mpi_rank == 0)
          std::cout << "Eigenvector orthonormality:\n";
        for (size_t root = 0; root < solver.n_roots(); root++) {
          for (size_t soot = 0; soot < solver.n_roots(); soot++) {
            auto result = handlers->rr().dot(x[root], x[soot]);
            if (mpi_rank == 0)
              std::cout << " " << result;
          }
          if (mpi_rank == 0)
            std::cout << std::endl;
        }
      }
      // if (mpi_rank == 0)
      //  std::cout << solver.statistics() << std::endl;
    }
  }
  MPI_Finalize();
}
