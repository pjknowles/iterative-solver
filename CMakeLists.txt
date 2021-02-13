cmake_minimum_required(VERSION 3.15)
project(LinearAlgebra CXX)

include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab+deploy-token-178345:RQZa2rmdwnV7Hdwg8XEq@gitlab.com/molpro/LibraryManager.git
        GIT_TAG 6d1737a0e35425c32cccec87741c1fbcc98751ac
)
FetchContent_MakeAvailable(LibraryManager)

LibraryManager_Project(FORTRAN_OPTION MPI_OPTION)

##############################
### Options

if (NOT MPI AND FORTRAN)
    message(FATAL_ERROR "Fortran serial implementation for LinearAlgebra is not yet ready") # TODO complete Fortran Serial implementation
endif ()
option(PPIDD "Whether to build on PPIDD" OFF)
option(LINEARALGEBRA_ARRAY_MPI3 "Build with MPI3 implementation of distributed array" ${MPI})
option(LINEARALGEBRA_ARRAY_GA "Build with GA implementation of distributed array" ${MPI})
option(LINEARALGEBRA_ARRAY_HDF5 "Build with HDF5 implementation of distributed array" ${MPI})
option(MPIFORTRAN "Build MPI fortran examples, using fortran MPI library" OFF)
if (PPIDD AND NOT PPIDD-GA)
    set(LINEARALGEBRA_ARRAY_GA OFF)
endif ()
##############################


add_subdirectory(dependencies)

add_subdirectory(src)

if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    add_subdirectory(benchmark)
    enable_testing()
    add_subdirectory(examples)
    add_subdirectory(test)
endif ()

add_subdirectory(doc)
