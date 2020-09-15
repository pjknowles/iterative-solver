add_executable(test-${PROJECT_NAME}-itsolv.exe testMatrix.cpp test_util.cpp)
target_link_libraries(test-${PROJECT_NAME}-itsolv.exe PUBLIC molpro::${PROJECT_NAME} gmock_main)
gtest_discover_tests(test-${PROJECT_NAME}-itsolv.exe)
