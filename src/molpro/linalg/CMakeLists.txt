# Extra checking for HDF5 and global arrays will go here
LibraryManager_Append(${PROJECT_NAME}
        SOURCES itsolv/IterativeSolver-double.cpp itsolv/IterativeSolver-complex-double.cpp options.cpp
        PUBLIC_HEADER IterativeSolverC.h itsolv/helper.h itsolv/ArrayHandlers.h itsolv/Statistics.h options.h
        PRIVATE_HEADER itsolv/helper-implementation.h
        )
if (FORTRAN)
    LibraryManager_Append(${PROJECT_NAME} SOURCES IterativeSolverF.F90 Iterative_Solver_Problem.F90 Iterative_Solver_PSpace.F90 Iterative_Solver_Matrix_Problem.F90)
endif ()

add_subdirectory(itsolv)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

if (MPI)
    message(STATUS "Building iterative-solver on MPI")
    LibraryManager_Append(${PROJECT_NAME} SOURCES IterativeSolverCMPI.cpp)
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAVE_MPI_H)
    target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)
else ()
    LibraryManager_Append(${PROJECT_NAME} SOURCES IterativeSolverCSerial.cpp)
    message(STATUS "Building serial iterative-solver")
endif ()

if (TARGET ppidd::ppidd)
    message(STATUS "Building iterative-solver on PPIDD library")
    target_link_libraries(${PROJECT_NAME} PUBLIC ppidd::ppidd)
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAVE_PPIDD_H)
endif ()

include(DependencyManager)
DependencyManager_Populate(profiler)
DependencyManager_Populate(utilities)
target_link_libraries(${PROJECT_NAME} PUBLIC molpro::profiler molpro::utilities)

LibraryManager_LAPACK(${PROJECT_NAME} PUBLIC Intel10_64lp_seq Intel10_64lp)
if (MKL_TYPE STREQUAL "lp64")
    message(STATUS "Eigen will use MKL")
    target_compile_definitions(${PROJECT_NAME} PUBLIC EIGEN_USE_MKL_ALL)
endif ()

if (NOT TARGET Eigen3::Eigen)
    message(VERBOSE "Eigen3::Eigen is not yet a target")
    find_package(Eigen3 3.3.0 QUIET)
    if (EIGEN3_FOUND)
        message(STATUS "Eigen3 found on system")
    else ()
        message(STATUS "Eigen3::Eigen not found on system, and will be downloaded")
        include(FetchContent)
        FetchContent_Declare(eigen3
                GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
                GIT_TAG 3.3.7
                )
        FetchContent_MakeAvailable(eigen3)
        add_library(Eigen3::Eigen ALIAS eigen)
    endif ()
else ()
    message(STATUS "iterative-solver: Eigen3::Eigen target already defined")
endif ()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    set(EIGEN_DONT_VECTORIZE ON)
endif ()
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake")
message("CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}")

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_MPL2_ONLY)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_DONT_VECTORIZE)
endif ()


include(LibraryManagerUtils)
print_target_properties(BLAS::BLAS STATUS)
print_target_properties(LAPACK::LAPACK )
print_target_properties(Eigen3::Eigen DEBUG)
print_target_properties(ga::ga DEBUG)
print_target_properties(ppidd::ppidd DEBUG)
print_target_properties(${PROJECT_NAME} DEBUG)

add_subdirectory(array)
