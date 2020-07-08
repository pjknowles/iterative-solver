cmake_minimum_required(VERSION 3.15)
project(LinearAlgebra CXX)

include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab+deploy-token-178345:RQZa2rmdwnV7Hdwg8XEq@gitlab.com/molpro/LibraryManager.git
        GIT_TAG master
)
FetchContent_MakeAvailable(LibraryManager)

LibraryManager_Project(FORTRAN_OPTION MPI_OPTION)

##############################
### Options

option(PPIDD "Whether to build on PPIDD" OFF)
option(LINEARALGEBRA_ARRAY_MPI3 "Build with MPI3 implementation of distributed array" ON)
option(LINEARALGEBRA_ARRAY_GA "Build with GA implementation of distributed array" ON)
##############################


add_subdirectory(dependencies)

add_subdirectory(src)

if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    add_subdirectory(examples)
    add_subdirectory(benchmark)
    enable_testing()
    add_subdirectory(test)
endif ()

# Documentation has to come last to make sure all dependencies have been added
include(DependencyManagerDocs)
DependencyManagerDocs_Add(${PROJECT_NAME}
        FILES README.md
        DOC_URL "https://molpro.gitlab.io/IterativeSolver/${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
        DEPENDS Profiler # projects whose documentation has to be built before current project.
        )
