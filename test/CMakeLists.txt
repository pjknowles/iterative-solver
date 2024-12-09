include(FetchContent)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "v1.14.0"
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

add_executable(test-${PROJECT_NAME}-test test.cpp)
target_link_libraries(test-${PROJECT_NAME}-test PUBLIC molpro::${PROJECT_NAME} gmock_main)
gtest_discover_tests(test-${PROJECT_NAME}-test)
if (VALGRIND)
  add_test(valgrind-${PROJECT_NAME}-test valgrind --error-exitcode=1 --leak-check=full test-${PROJECT_TEST}-test)
endif ()



add_subdirectory(array)
add_subdirectory(itsolv)
