iterative-solver
================

[![Build and test](https://github.com/molpro/iterative-solver/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/molpro/iterative-solver/actions/workflows/build-and-test.yml)

[//]: # (&#40;https://github.com/molpro/iterative-solver/commits/master&#41;)
[![license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/molpro/iterative-solver/blob/master/LICENSE)

[![Online documentation - Docs](https://img.shields.io/badge/Online_documentation-Docs-2ea44f)](https://molpro.github.io/iterative-solver/)

## Overview

Implements iterative solvers for linear and non-linear problems and
distributed arrays for HPC. The solvers are specialised to work with specific data types, but are also templated on the
container allowing for easy integration into existing software.

List of key features:
* Implements iterative solvers for eigenvalue problem, linear equations, optimisation (L-BFGS) and non-linear equations 
(DIIS)
* Novel algorithms including P-space and D-space (see paper)
* Structured to allow easy addition of new solvers or modification of the current ones without changes to the user's 
code
* Templated on container for easy integration into existing programs
  *  User defined containers can be used without modification with the help of array handler abstraction
* Precision-agnostic: the dense kernels dispatch to BLAS/LAPACK for the scalar types LAPACK covers and to
Eigen for any other, so the solvers can be used in extended or arbitrary precision
* Explicitly instantiated for double, complex double and long double value types, so that all heavy numerical
operations are only compiled once for those
* Provides distributed arrays in memory and on disk for HPC
* Contains Fortran and C wrappers

## Installation

CMake is used to build the library and it can integrate easily with other CMake builds.

```cmake
include(FetchContent)
FetchContent_Declare(
        iterative-solver
        GIT_REPOSITORY https://github.com/molpro/iterative-solver.git
        GIT_TAG ${COMMIT_HASH_OR_TAG_VALUE})
FetchContent_MakeAvailable(iterative-solver)
target_link_libraries(${YOUR_LIBRARY_NAME} PUBLIC molpro::iterative-solver)
```

## Usage

### Interfaces

There is a hierarchy of abstract classes defined in `molpro/linalg/itsolv/IterativeSolver.h` with `IterativeSolver` 
defining the interface for common functionality and `LinearEigensystem`, `LinearEquations`, `Optimize` and
`NonLinearEquations` defining functions that are specific to each type of solver. They are provided to reduce header 
bloat in user's code. 

Each type of solver has at least one implementation class, e.g. `LinearEigensystemDavidson` in 
`molpro/linalg/itsolv/LinearEigensystemDavidson.h`, which implements the full interface. The library is designed to make it easy
to add new iterative solvers, so there might be more than one implementation.

### Example

The simplest way to use the library is to define a `Problem` subclass that supplies
the matrix action and (optionally) diagonals, then call `solve()`. Here is an example
of using `LinearEigensystemDavidson`,

```cpp
#include <molpro/linalg/itsolv/LinearEigensystemDavidson.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
// ...
using R = std::vector<double>;
class MyProblem : public molpro::linalg::itsolv::Problem<R> {
    void action(const CVecRef<R>& parameters, const VecRef<R>& actions) const override { /* ... */ }
    bool diagonals(R& d) const override { /* fill d */ return true; }
};

auto solver = molpro::linalg::itsolv::create_LinearEigensystem<R>();
R parameters(n), actions(n);
MyProblem problem;
if (!solver->solve(parameters, actions, problem, /*generate_initial_guess=*/true)) {
    // deal with the unconverged case
}
solver->solution(parameters, actions);
```

More advanced users can copy and modify the code in `solve()` to tailor it for their own use, see implementation in `molpro/linalg/itsolv/IterativeSolverTemplate.h`.

### Containers and array handlers

In many programs there are special containers for storing the vectors and operating on them. This might be for efficiency reasons,
e.g. exploiting symmetry of the problem, or for collecting metadata, e.g. memory usage and operation count. In either case,
The code is templated on container types to ease adaptation.

To avoid the code-bloat of header only libraries all of the numerically intensive work is explicitly instantiated for
`double`, `std::complex<double>` and `long double` element types.
This makes recompilation of IterativeSolver with different container types very fast.
Containers with other element types are supported as well; see [Arbitrary precision](#arbitrary-precision).

There are 3 types of containers:

* R &mdash; working set container
  * fast access to elements
  * used sparingly to preserve memory 
* Q &mdash; slow container
  * slow access to elements
  * used for most of the subspace
  * usually on disk where there is lots of storage space
* P &mdash; sparse container
  * light-weight sparse container (e.g. `std::map<size_t, double>`)

IterativeSolver does not modify containers directly instead using ArrayHandler for all array related operations. This allows for only modest
restrictions on containers:

* must have a  move constructor

* must have an element type the solvers can compute in; see [Arbitrary precision](#arbitrary-precision)

ArrayHandler is an abstract class used by IterativeSolver to perform copy and linear algebra operations (`dot`, `axpy`). 
iterative-solver provides implementations for Iterable containers (e.g. `std::vector`), distributed containers (e.g. `molpro::linalg::array::DistrArray`),
and mapped containers (e.g. `std::map`). However, some users might need/want to provide their own implementations. 

### Arbitrary precision

Nothing in the solvers is tied to double precision. The scalar type is taken from the container
(`R::value_type`), and every threshold, working array and dense decomposition follows it.

The dense linear algebra dispatches on that scalar type: `float`, `double`, `std::complex<float>` and
`std::complex<double>` go to LAPACK through the LAPACKE interface, and every other type goes to the
equivalent templated Eigen decomposition. That covers any scalar for which `Eigen::NumTraits` is
specialised — `long double` out of the box, and extended- and arbitrary-precision types such as
`boost::multiprecision` or `mpfr::mpreal` through the Eigen support those libraries provide. The
dispatch is in `molpro/linalg/itsolv/helper-dispatch.h`; `has_lapack_kernel<T>` reports which branch a
given type takes.

```cpp
using R = std::vector<long double>;
auto solver = molpro::linalg::itsolv::create_LinearEigensystem<R>();
solver->set_convergence_threshold(1e-17L); // out of reach in double precision
```

The library ships explicit instantiations of the dense kernels for `long double` alongside those for
`double` and `std::complex<double>`, so that case costs no extra compilation. For any other scalar
type the kernels are instantiated implicitly from `molpro/linalg/itsolv/helper-implementation.h`. The
solver classes themselves are instantiated for your container types as usual, by writing
`template class SolverFactory<R, Q, P>;` in one translation unit that includes
`molpro/linalg/itsolv/SolverFactory-implementation.h`.

Thresholds that were originally hard-coded for double precision — "this quantity is zero", "this
direction is null" — are not absolute constants but are rescaled to the working precision by
`precision_scaled<T>()`, which preserves their margin measured in machine epsilons and returns them
unchanged when the working precision is double. Without this the solvers would declare convergence,
or discard new directions as redundant, long before the extra precision had been used.

One limitation remains: the distributed arrays in `molpro::linalg::array` fix `value_type` to
`double`, so extended precision applies to in-memory containers such as `std::vector<long double>`.

### Complex arithmetic

Containers with a complex element type are supported throughout, and compose with the above: a
`std::vector<std::complex<long double>>` takes the Eigen branch of the dispatch at extended precision.

**The inner product is hermitian.** `ArrayHandler::dot(x, y)` is `<x|y>` -- conjugate-linear in its
first argument, linear in its second -- so `dot(x, x)` is real and non-negative and the subspace
overlap is a hermitian matrix. Earlier releases computed `sum(x_i * y_i)` without the conjugation,
which made the overlap complex symmetric; that was never usable, because every complex subspace solver
was an `assert(false)` stub, but the change matters to anyone who wrote their own `ArrayHandler`
subclass for a complex container. If you want the bilinear c-product of complex-scaling and CAP
methods instead, do not conjugate in your own handler -- but note that the subspace solvers assume the
hermitian convention.

Quantities that are magnitudes rather than scalars of the problem -- residual norms, singular values,
convergence thresholds -- are typed `value_type_abs`, which is `double` for a `std::complex<double>`
container. `IterativeSolver::errors()` returns them, and so does `set_convergence_threshold()`.

The optimisers (`Optimize`, i.e. BFGS and steepest descent) accept complex containers, but their
objective is taken to be real: the line search orders function values, so only the real part of the
value returned by `Problem::residual()` is used. Selecting a P space also needs an ordering and so
throws for a complex element type, as it did before.

## Contributing

[doc/contributing/ConsistencyAudit.md](doc/contributing/ConsistencyAudit.md) records the playbook
used for the consistency-audit series, for anyone wanting to redo or extend that work.

## Citing

Any publications resulting from this work should cite relevant papers in CITE.txt

## List of Contributors

Peter Knowles

Marat Sibaev

Iakov Polyak

Rob Welch
