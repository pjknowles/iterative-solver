include(FetchContent)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.10.0"
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)
foreach (t IterativeSolver-test IterativeSolver-OOCA-test TestIterativeSolver PagedVector-test OutOfCoreArray-test SimpleVector-test)
    add_executable(${t}.exe ${t}.cpp test.h)
    target_link_libraries(${t}.exe PUBLIC molpro::${PROJECT_NAME} gmock_main)# ${MPI_CXX_LIBRARIES})
    gtest_discover_tests(${t}.exe)
    message(STATUS "configuring test ${t}")
    include(LibraryManagerUtils)
    print_variable(MPIEXEC_EXECUTABLE)
    print_variable(MPIPROCS)
    if (MPI_FOUND)
        if (NOT MPIPROCS)
            set(MPIPROCS 1)
        endif ()
        add_test(${t} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIPROCS} ${MPIEXEC_PREFLAGS} ${t}.exe ${MPIEXEC_POSTFLAGS})
    else (MPI_FOUND)
        add_test(${t} ${t}.exe)
    endif (MPI_FOUND)
endforeach ()

add_subdirectory(array)
