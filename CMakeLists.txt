cmake_minimum_required(VERSION 3.15)
project(LinearAlgebra CXX)

include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab+deploy-token-178345:RQZa2rmdwnV7Hdwg8XEq@gitlab.com/molpro/LibraryManager.git
        GIT_TAG fef21c35ac00abea4d810ee8c8dbc34dcf8fa2e8
)
FetchContent_MakeAvailable(LibraryManager)

LibraryManager_Project(FORTRAN_OPTION MPI_OPTION)

##############################
### Options

option(PPIDD "Whether to build on PPIDD" OFF)
option(LINEARALGEBRA_ARRAY_MPI3 "Build with MPI3 implementation of distributed array" ON)
option(LINEARALGEBRA_ARRAY_GA "Build with GA implementation of distributed array" ON)
option(LINEARALGEBRA_ARRAY_HDF5 "Build with HDF5 implementation of distributed array" ON)
if (PPIDD)
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
