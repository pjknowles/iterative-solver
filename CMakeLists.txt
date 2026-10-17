cmake_minimum_required(VERSION 3.5...4.0)
project(iterative-solver CXX)

# The C++20 requirement makes CMake scan every source for modules, which needs
# clang-scan-deps alongside clang, and Debian's clang package does not provide it.
# The project has no modules, so the scan only costs time or fails the build.
if (NOT DEFINED CMAKE_CXX_SCAN_FOR_MODULES)
    set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
endif ()

include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://github.com/molpro/library-manager.git
        GIT_TAG 0.13.5
)
FetchContent_MakeAvailable(LibraryManager)

LibraryManager_Project(FORTRAN_OPTION MPI_OPTION DEFAULT_VERSION "$Format:%(describe:tags)$")

##############################
### Options

if (NOT MPI AND FORTRAN)
    message(FATAL_ERROR "Fortran serial implementation for iterative-solver is not yet ready") # TODO complete Fortran Serial implementation
endif ()
option(PPIDD "Whether to build on PPIDD" OFF)
option(LINEARALGEBRA_ARRAY_MPI3 "Build with MPI3 implementation of distributed array" ${MPI})
option(LINEARALGEBRA_ARRAY_GA "Build with GA implementation of distributed array" ${MPI})
option(LINEARALGEBRA_ARRAY_HDF5 "Build with HDF5 implementation of distributed array" OFF)
option(MPIFORTRAN "Build MPI fortran examples, using fortran MPI library" OFF)
option(LIBRARY_ONLY "Build only the library" OFF)
option(BENCHMARK "Build benchmark executables" OFF)
option(VALGRIND "Run valgrind tests and examples" OFF)
find_program(VALGRIND_EXISTS "valgrind")
if (NOT VALGRIND_EXISTS)
    set(VALGRIND OFF)
endif ()
if (PPIDD AND NOT PPIDD-GA)
    set(LINEARALGEBRA_ARRAY_GA OFF)
endif ()
message(STATUS "iterative-solver: MPI=${MPI}")
message(STATUS "iterative-solver: FORTRAN=${FORTRAN}")
message(STATUS "iterative-solver: LINEARALGEBRA_ARRAY_MPI3=${LINEARALGEBRA_ARRAY_MPI3}")
message(STATUS "iterative-solver: LINEARALGEBRA_ARRAY_GA=${LINEARALGEBRA_ARRAY_GA}")
message(STATUS "iterative-solver: LINEARALGEBRA_ARRAY_HDF5=${LINEARALGEBRA_ARRAY_HDF5}")
##############################

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    # CMAKE_<LANG>_COMPILER_LAUNCHER alone is sufficient for both the Makefiles and
    # Ninja generators. Also setting RULE_LAUNCH_COMPILE double-wraps every compile
    # command (ccache ccache <compiler> ...), which confuses ccache's own argument
    # parsing and can produce spurious "file not found" errors.
    set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()


add_subdirectory(dependencies)

add_subdirectory(src)

if (NOT LIBRARY_ONLY)
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    # benchmark.cpp calls molpro::mpi::init()/finalize(), and molpro/mpi.h keys HAVE_MPI_H off
    # __has_include(<mpi.h>) rather than this project's MPI option, so in a serial build on a
    # machine that has MPI headers it compiles MPI calls that are never linked.
    if (MPI)
        add_subdirectory(benchmark)
    endif ()
    enable_testing()
    add_subdirectory(examples)
    add_subdirectory(test)
endif ()

add_subdirectory(doc)
endif()
