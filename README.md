Welcome to LinearAlgebra
========================

[![pipeline status](https://gitlab.com/molpro/linearalgebra/badges/master/pipeline.svg)](https://gitlab.com/molpro/linearalgebra/-/commits/master)
[![license](https://img.shields.io/badge/license-MIT-green.svg)](https://gitlab.com/molpro/linearalgebra/-/blob/master/LICENSE)
[![license](https://img.shields.io/badge/documentation-blue.svg)](https://dependencymanager.gitlab.io/dependency-manager/)


## Overview

LinearAlgebra is a part of Molpro's ecosystem. It implements iterative solvers for linear and non-linear problems and
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
* Specialised for double and complex value types, so that all heavy numerical operations are only compiled once
* Provides distributed arrays in memory and on disk for HPC
* Contains Fortran and C wrappers

## Installation

CMake is used to build the library and it can integrate easily with other CMake builds.

```cmake
include(FetchContent)
FetchContent_Declare(
        linearalgebra
        GIT_REPOSITORY https://gitlab.com/molpro/linearalgebra.git
        GIT_TAG ${COMMIT_HASH_OR_TAG_VALUE})
FetchContent_MakeAvailable(linearalgebra)
target_link_libraries(${YOUR_LIBRARY_NAME} PUBLIC molpro::LinearAlgebra)
```

## Usage

### Interfaces

There is a hierarchy of abstract classes defined in `molpro/linalg/istolv/IterativeSolver.h` with `IterativeSolver` 
defining the interface for common functionality and `ILinearEigensystem`, `ILinearEquations`, `IOptimize` and
`INonLinearEquations` defining functions that are specific to each type of solver. They are provided to reduce header 
bloat in user's code. 

Each type of solver has at least one implementation class, e.g. `LinearEigensystem` in 
`molpro/linalg/itsolv/LinearEigensystem.h`, which implements the full interface. The library is designed to make it easy
to add new iterative solvers, so there might be more than one implementation.

### Example

The simplest way to use the library is to call `solve()` and pass a function for getting the action of the matrix on parameter set. 
Here is an example of using `LinearEigensystem` with Davidson preconditioner,

```cpp
#include <molpro/linalg/itsolv/LinearEigensystem.h>
// ...
using molpro::linalg::itsolv::LinearEigensystem;
using R = std::vector<double>;
using Q = std::vector<double>;
LinearEigensystem<R, Q> solver{};
auto [x, g] = get_initial_parameters_and_action();
solver.set_preconditioner_davidson(get_diagonal_matrix_elements());
solver.solve(x, g, user_specified_action_function);
if (!solver.converged()){
    // deal with the unconverged case
}
```

More advanced users can copy and modify the code in `solve()` to tailor it for their own use, see implementation in `molpro/linalg/itsolv/IterativeSolverTemplate.h`.

### Containers and array handlers

In many programs there are special containers for storing the vectors and operating on them. This might be for efficiency reasons,
e.g. exploiting symmetry of the problem, or for collecting metadata, e.g. memory usage and operation count. In either case,
IterativeSolver is templated on container types to ease adaptation.

To avoid the code-bloat of header only libraries all of the numerically intensive work is specialised for `double` and `std::complex<double>` types.
This makes recompilation of IterativeSolver with different container types very fast.

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

* must have conforming element type (currently `double` and `std::complex<double>`, but this can be easily extended)

ArrayHandler is an abstract class used by IterativeSolver to perform copy and linear algebra operations (`dot`, `axpy`). 
LinearAlgebra provides implementations for Iterable containers (e.g. `std::vector`), distributed containers (e.g. `molpro::linalg::array::DistrArray`),
and mapped containers (e.g. `std::map`). However, some users might need/want to provide their own implementations. 

## Citing

Any publications resulting from this work should cite relevant papers in CITE.txt

## Contributing

This library is under active development and potential collaborators are welcome to contact
Prof. Peter Knowles at Cardiff University.

## List of Contributors

Prof. Peter Knowles

Marat Sibaev

Iakov Polyak
