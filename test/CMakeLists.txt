include(FetchContent)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.10.0"
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)
foreach (t IterativeSolver-test TestIterativeSolver PagedVector-test OutOfCoreArray-test SimpleVector-test LinearEigensystem-file-test)
    add_executable(${t}.exe ${t}.cpp test.h)
    target_link_libraries(${t}.exe PUBLIC molpro::${PROJECT_NAME} gmock_main)# ${MPI_CXX_LIBRARIES})
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
file(COPY ../examples/bh.hamiltonian DESTINATION .)

add_subdirectory(array)
add_subdirectory(itsolv)
