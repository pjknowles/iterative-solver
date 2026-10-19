# AGENTS.md

Guidance for AI agents (and humans) working on the iterative-solver repository.

## Code layout

- `src/molpro/linalg/` – the library (`molpro::iterative-solver`)
  - `itsolv/` – iterative solvers (namespace `molpro::linalg::itsolv`)
    - `IterativeSolver.h` – abstract interfaces (`IterativeSolver`, `LinearEigensystem`, `LinearEquations`,
      `Optimize`, `NonLinearEquations`) and `Problem`.
    - `IterativeSolverTemplate.h` – shared implementation of the interfaces, including the `solve()` driver.
    - One header per solver implementation (`LinearEigensystemDavidson.h`, `OptimizeBFGS.h`, ...), each paired with
      an options class (`<Solver>Options.{h,cpp}`) whose hierarchy mirrors the solver hierarchy.
    - `SolverFactory.h` / `SolverFactory-implementation.h` – dynamic creation of solvers from options.
    - `subspace/` – P/Q/R/D/X subspace management and subspace solvers.
    - `helper.h`, `helper-implementation.h`, `helper-dispatch.h` – dense linear algebra kernels (BLAS/LAPACK or Eigen).
    - `IterativeSolver-{double,complex-double,long-double}.cpp` – explicit instantiations of the dense kernels.
  - `array/` – array handlers (`ArrayHandler*.h`) and distributed arrays (`DistrArray*`: MPI3, GA, HDF5, file,
    span) in namespace `molpro::linalg::array`.
  - `IterativeSolverC.h`, `IterativeSolverCMPI.cpp`, `*.F90` – C and Fortran interfaces.
- `test/` – GoogleTest unit tests: `test/array/`, `test/itsolv/`, `test/itsolv/subspace/`.
- `examples/` – C++ and Fortran example programs (also run as tests) and Hamiltonian data files.
- `benchmark/` – array-operation benchmark (`benchmark.exe`), see [Benchmarks](#benchmarks).
- `python/` – Cython bindings (`iterative_solver`) and their `unittest` tests in `python/test/`.
- `doc/` – Doxygen sources (`*.dox`), `LibraryDesign.md` and `contributing/` (e.g. the consistency-audit playbook).
- `dependencies/` – fetching of external dependencies via the molpro LibraryManager/DependencyManager.

## Building and testing

- CMake build, C++20. Relevant options: `MPI`, `FORTRAN`, `FORTRAN_INTEGER8`, `LINEARALGEBRA_ARRAY_MPI3`,
  `LINEARALGEBRA_ARRAY_GA`, `LINEARALGEBRA_ARRAY_HDF5`, `BENCHMARK`, `LIBRARY_ONLY`, `VALGRIND`.
  `FORTRAN=ON` requires `MPI=ON`; a serial build is `-DMPI=OFF -DFORTRAN=OFF`.
- Always configure with `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON`.
- Tests are registered with `ctest` (`gtest_discover_tests` plus the examples). Run with
  `ctest --output-on-failure` from the build directory; use `-R <regex>` to select subsets.
- With MPI, solver tests run with `$MPIPROCS` processes (default 1); Fortran-backed tests additionally run as
  `*_3procs` (and some `*_2procs`) to exercise uneven distributions.
- Platform-dependent tests are skipped with `GTEST_SKIP()` and a reason (e.g. when `long double` is not wider than
  `double`), never silently passed.
- Python bindings: `python/build.sh`, then `python -m unittest discover -s python/test` (see `python/test-linux.sh`).
- CI (`.github/workflows/build-and-test.yml`) builds and tests a matrix of gcc/clang, MPICH/OpenMPI, serial, HDF5,
  OpenBLAS/MKL, macOS arm64, Linux aarch64 and an ASan/UBSan build. Changes should be expected to pass on all of
  them; in particular, do not rely on `long double` being 80-bit, or on `MPI_Comm` being an integer or a pointer.

## Benchmarks

- `benchmark/benchmark.exe` is built with every non-`LIBRARY_ONLY` build of this project, independent of the
  `BENCHMARK` option, and is registered as the ctest test `benchmark` (run with `$MPIPROCS` processes under MPI).
  Exclude it with `ctest -E '^benchmark$'` when only correctness matters, and run it on an otherwise idle machine
  when timings matter.
- It measures the array-handler operations (`dot`, `axpy`, `scal`, `fill`, `copy`, inner/outer `gemm`, and their
  transposed / both-direction variants) between a fast in-memory container (`DistrArraySpan`, `std::vector`, `Span`)
  and a slow one (`DistrArrayFile`, and `DistrArrayMPI3` / `DistrArrayHDF5` when enabled) for vector lengths from
  500 to 10^6, plus generic STL-iterator access versus native access of distributed arrays.
- The benchmark classes live in `benchmark/ArrayBenchmark.h` (`ArrayBenchmark` and factory functions such as
  `ArrayBenchmarkDDisk`, `ArrayBenchmarkIterableDDisk`). Each operation is repeated until it runs for roughly the
  given target time; results are reported via `molpro::Profiler` and printed to stdout, and some runs also write
  profiler call graphs (`*.gv`) into the working directory.
- `-DBENCHMARK=ON` together with `FORTRAN` additionally builds `LinearEigensystemBenchmarkF` (from
  `examples/LinearEigensystemBenchmarkF.F90`), an end-to-end Davidson benchmark through the Fortran interface. It is
  not a ctest test and has to be run by hand.
- When a change is meant to improve performance, compare benchmark output before and after the change with the same
  build type (Release) and the same number of MPI processes. Ensure the affected code is actually meaningfully
  exercised by the relevant benchmark(s) beforehand. If not, extend the benchmark set.

## Development guidelines

### Priorities

1. Correctness.
2. Computational efficiency, closely followed by readability and maintainability. If the code can be made more
   performant by writing it differently, do so, but weigh cost against benefit: a large performance gain can justify a
   large loss of readability, a negligible gain does not justify any.

### Design

- Avoid code duplication. Whenever possible, factor common functionality out into a shared helper function or class
  and reuse it at every call site, so that each piece of functionality has a single source of truth.
- The code may run with different value types: `double`, `std::complex<double>`, `long double` and, via the Eigen
  fallback, other (e.g. arbitrary-precision) scalar types. Any code reachable with one of these types must work
  correctly for it. For complex values in particular, mind conjugation in inner products, use `real_type_t<T>` for
  norms, thresholds and other real-valued quantities, and do not compare or order complex values directly.
- Numerically heavy, value-type-dependent functions are declared `extern template` in the headers and explicitly
  instantiated for `double`, `std::complex<double>` and `long double` in `IterativeSolver-*.cpp`. A new kernel needs
  all three.
- Solvers are templated on the containers `R`, `Q`, `P`; all container operations go through `ArrayHandler`s, never
  directly on the container.
- New solver options go into the corresponding `*Options` class and must be wired into `set_options()` /
  `get_options()` (and `SolverFactory` for new solvers).
- Errors: invalid arguments in array handlers are reported via the handler's `error()` (`ArrayHandlerError`);
  elsewhere throw standard exceptions with a descriptive message. Check preconditions before modifying output
  arguments. Destructors must not throw.
- Diagnostic output goes through the solver's `Logger`, not `std::cout`; verbosity-dependent output is gated by the
  solver verbosity.
- Keep the profiler balanced: prefer the RAII handle returned by `profiler()->push(...)`, and make sure every
  `start()` has a matching `stop()`, including on early returns.

### Style and documentation

- Format with the repository's `.clang-format` (LLVM based, 120 columns, left-aligned pointers/references).
- Include guards (not `#pragma once`) named after the path, e.g.
  `LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMDAVIDSON_H`.
- Includes of library headers use angle brackets and the full path: `#include <molpro/linalg/itsolv/helper.h>`.
- Nested namespace definitions (`namespace molpro::linalg::itsolv {`), closed with `} // namespace ...`.
  No file-scope `using` declarations in headers.
- Keep in-source comments brief and to the point, and describe only the code as it currently is. Larger explanations,
  the reasoning behind a change, and comparisons with what used to be there belong in the commit message. The only
  exception: code that was deliberately changed into a non-obvious form may carry a brief reason, so that nobody
  "fixes" it back to the old version.
- Document public API with Doxygen comments (`/*! @brief ... */`, `@param`, `@copydoc` for overloads that behave the
  same). Keep the `.dox` files in `doc/` and the README in sync with behaviour changes.

### Tests

- Every bug fix and feature comes with a GoogleTest test in the matching `test/` subdirectory; add new test files to
  the target's source list in the corresponding `CMakeLists.txt`.
- Code reachable with several value types is tested for them; extend the type-specific tests (e.g.
  `test_helper-complex.cpp`, `test_complex_e2e.cpp`, `test_precision_e2e.cpp`) accordingly.
- Distributed-array and Fortran code must also be correct with more than one MPI process.

### Commits and pull requests

- Each commit is one logical changeset: all changes that logically belong together go into the same commit, and
  changes that are logically different go into separate commits.
- Subject line in imperative mood, no trailing period, describing the effect
  (e.g. "Check sparse indices before overwriting the destination of copy"). CI-only changes are prefixed `CI: `.
- The commit body explains *why*: the faulty behaviour, how it arises (referencing commits by short hash and
  issues/PRs by `#number`), what was changed and which alternatives were rejected.
- Commits that were produced with AI assistance end with the line `AI-assisted` - do not include tool-specific
  attribution.
- Work happens on feature branches off `master`, merged via pull request.
