add_executable(test-${PROJECT_NAME}-itsolv.exe test_wrap.cpp test_util.cpp testDSpaceResetter.cpp test_SolverFactory.cpp
        test_options_map.cpp test_Options.cpp SolverFactory_implementation.cpp
        test_Interpolate.cpp
        test_RSPT.cpp
        )
target_link_libraries(test-${PROJECT_NAME}-itsolv.exe PUBLIC molpro::${PROJECT_NAME} gmock_main)
gtest_discover_tests(test-${PROJECT_NAME}-itsolv.exe)

file(COPY ../../examples/bh.hamiltonian DESTINATION .)
file(COPY ../../examples/hf.hamiltonian DESTINATION .)
file(COPY ../../examples/phenol.hamiltonian DESTINATION .)
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
    message(STATUS "LinearAlgebra: Eigen3::Eigen target already defined")
endif ()

foreach (t test_LinearEigensystem test_LinearEquations test_Optimize)
    add_executable(${t}.exe ${t}.cpp test.h SolverFactory_implementation.cpp)
    if (FORTRAN)
        target_sources(${t}.exe PRIVATE ${t}F.f90)
    else()
        target_compile_definitions(${t}.exe PRIVATE -DNOFORTRAN)
    endif()
    target_link_libraries(${t}.exe PUBLIC molpro::${PROJECT_NAME} gmock_main)# ${MPI_CXX_LIBRARIES})
    target_link_libraries(${t}.exe PUBLIC Eigen3::Eigen)
    gtest_discover_tests(${t}.exe)
    if (MPI)
        # in principle this should be configured at runtime but it seems cmake/ctest support is not yet there to do it simply
        if (NOT DEFINED ENV{MPIPROCS})
            set(ENV{MPIPROCS} 1)
        endif ()
        message(DEBUG "add_test(${t} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} $ENV{MPIPROCS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${t}.exe ${MPIEXEC_POSTFLAGS})")
        add_test(${t} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} $ENV{MPIPROCS} ${MPIEXEC_PREFLAGS} "${CMAKE_CURRENT_BINARY_DIR}/${t}.exe" ${MPIEXEC_POSTFLAGS})
    else ()
        add_test(${t} ${t}.exe)
    endif ()
endforeach ()

add_subdirectory(subspace)