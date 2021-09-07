include(FetchContent)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.11.0"
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

add_executable(test-${PROJECT_NAME}-test test.cpp)
target_link_libraries(test-${PROJECT_NAME}-test PUBLIC molpro::${PROJECT_NAME} gmock_main)
gtest_discover_tests(test-${PROJECT_NAME}-test)



add_subdirectory(array)
add_subdirectory(itsolv)
