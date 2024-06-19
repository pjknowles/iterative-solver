cmake_minimum_required(VERSION 3.15)
project(iterative-solver CXX)

include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://github.com/molpro/library-manager.git
        GIT_TAG 0.13.0
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
if (PPIDD AND NOT PPIDD-GA)
    set(LINEARALGEBRA_ARRAY_GA OFF)
endif ()
message(STATUS "iterative-solver: MPI=${MPI}")
message(STATUS "iterative-solver: FORTRAN=${FORTRAN}")
message(STATUS "iterative-solver: LINEARALGEBRA_ARRAY_MPI3=${LINEARALGEBRA_ARRAY_MPI3}")
message(STATUS "iterative-solver: LINEARALGEBRA_ARRAY_GA=${LINEARALGEBRA_ARRAY_GA}")
message(STATUS "iterative-solver: LINEARALGEBRA_ARRAY_HDF5=${LINEARALGEBRA_ARRAY_HDF5}")
##############################


add_subdirectory(dependencies)

add_subdirectory(src)

if (NOT LIBRARY_ONLY)
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    add_subdirectory(benchmark)
    enable_testing()
    add_subdirectory(examples)
    add_subdirectory(test)
endif ()

add_subdirectory(doc)
endif()
