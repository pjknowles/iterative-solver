foreach (t IterativeSolver-test PagedVector-test SimpleVector-test)
    add_executable(${t} ${t}.cpp)
    target_link_libraries(${t} ${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
    set_target_properties(${t} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
    if (MPI_FOUND)
        add_test(${t} ${MPIEXEC} ${MPIOPTIONS} ${PROJECT_BINARY_DIR}/${t})
    else (MPI_FOUND)
        add_test(${t} ${PROJECT_BINARY_DIR}/${t})
    endif (MPI_FOUND)
endforeach ()

include(FetchContent)
FetchContent_Declare(
        googletest_dep
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.8.1"
)
FetchContent_GetProperties(googletest_dep)
if (NOT googletest_dep_POPULATED)
    FetchContent_Populate(googletest_dep)
    add_subdirectory(${googletest_dep_SOURCE_DIR} ${googletest_dep_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()


add_executable(test.exe TestIterativeSolver.cpp)
target_link_libraries(test.exe PUBLIC ${PROJECT_NAME} gmock_main)

include(GoogleTest)
gtest_discover_tests(test.exe)