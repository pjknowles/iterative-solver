if (MPI)
    add_executable(test-${PROJECT_NAME}-parallel
            main.cpp
            testutil.cpp
            testDistrArray.h
            testDistrArrayMPI3.cpp
            testDistrArrayGA.cpp
            testHDF5Handle.cpp
            )
    target_link_libraries(test-${PROJECT_NAME}-parallel PUBLIC molpro::${PROJECT_NAME} gmock_main)
    target_compile_definitions(test-${PROJECT_NAME}-parallel PUBLIC "ARRAY_DATA=\"${CMAKE_CURRENT_SOURCE_DIR}/data\"")
    gtest_discover_tests(test-${PROJECT_NAME}-parallel)

    add_test(NAME test-${PROJECT_NAME}-parallel
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIOPTIONS} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            ${CMAKE_CURRENT_BINARY_DIR}/test-${PROJECT_NAME}-parallel ${MPIEXEC_POSTFLAGS}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            )
endif ()
