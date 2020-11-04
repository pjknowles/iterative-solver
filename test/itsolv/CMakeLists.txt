add_executable(test-${PROJECT_NAME}-itsolv.exe test_wrap.cpp test_util.cpp test_reset_dspace.cpp)
target_link_libraries(test-${PROJECT_NAME}-itsolv.exe PUBLIC molpro::${PROJECT_NAME} gmock_main)
gtest_discover_tests(test-${PROJECT_NAME}-itsolv.exe)

add_subdirectory(subspace)