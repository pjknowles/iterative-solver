LibraryManager_Append(${PROJECT_NAME}
        SOURCES
        IterativeSolver.cpp
        PUBLIC_HEADER
        SimpleVector.h
        PagedVector.h
        OpaqueVector.h
        OutOfCoreArray.h
        IterativeSolver.h
        )
if (CMAKE_Fortran_COMPILER)
    LibraryManager_Append(${PROJECT_NAME}
            SOURCES IterativeSolverF.F90
            )
endif ()

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

include(DependencyManager)
DependencyManager_Populate(Profiler)
DependencyManager_Populate(utilities)
target_link_libraries(${PROJECT_NAME} PUBLIC molpro::Profiler molpro::utilities)

if (MPI_CXX_FOUND)
    message(STATUS "Building IterativeSolver on MPI")
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAVE_MPI_H)
    target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)
else ()
    message(STATUS "Building serial IterativeSolver")
endif ()
if (TARGET ppidd::ppidd)
    message(STATUS "Building IterativeSolver on PPIDD library")
    target_link_libraries(${PROJECT_NAME} PUBLIC ppidd::ppidd)
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAVE_PPIDD_H)
endif ()

LibraryManager_LAPACK(${PROJECT_NAME} PUBLIC Intel10_64lp_seq Intel10_64lp)
if (MKL_TYPE STREQUAL "lp64")
    message(STATUS "Eigen will use MKL")
    target_compile_definitions(${PROJECT_NAME} PUBLIC EIGEN_USE_MKL_ALL)
endif ()

if (NOT TARGET Eigen3::Eigen)
    message(VERBOSE "Eigen3::Eigen is not yet a target")
    find_package(Eigen3 3.3.0)
    if (EIGEN3_FOUND)
        message(STATUS "Eigen3 found on system")
    else ()
        if (NOT CMAKE_Fortran_COMPILER)
            message(FATAL_ERROR "Eigen3 cannot be configured from source with FORTRAN=OFF; install it in the system first")
        endif ()
        include(FetchContent)
        FetchContent_Declare(eigen3
                GIT_REPOSITORY https://gitlab.com/libeigen/eigen
                GIT_TAG 3.3.7
                )
        FetchContent_MakeAvailable(eigen3)
        message(STATUS "Eigen3::Eigen not found on system, and downloaded")
    endif ()
else ()
    message(STATUS "SymmetryMatrix: Eigen3::Eigen target already defined")
endif ()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    set(EIGEN_DONT_VECTORIZE ON)
endif ()

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_DONT_VECTORIZE)
endif ()


include(LibraryManagerUtils)
print_target_properties(LAPACK::LAPACK DEBUG)
print_target_properties(Eigen3::Eigen DEBUG)
print_target_properties(${PROJECT_NAME} DEBUG)