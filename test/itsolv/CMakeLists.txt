add_executable(test-${PROJECT_NAME}-itsolv.exe test_wrap.cpp test_util.cpp testDSpaceResetter.cpp test_SolverFactory.cpp
        test_options_map.cpp test_Options.cpp SolverFactory_implementation.cpp Optimize-test.cpp)
target_link_libraries(test-${PROJECT_NAME}-itsolv.exe PUBLIC molpro::${PROJECT_NAME} gmock_main)
gtest_discover_tests(test-${PROJECT_NAME}-itsolv.exe)
file(COPY ../../examples/bh.hamiltonian DESTINATION .)
file(COPY ../../examples/hf.hamiltonian DESTINATION .)
file(COPY ../../examples/phenol.hamiltonian DESTINATION .)
add_subdirectory(subspace)