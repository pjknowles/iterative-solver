if (MPI)
    if (LINEARALGEBRA_ARRAY_HDF5)
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
        foreach (fil empty_file inner_group_dataset.hdf5 single_dataset.hdf5)
            configure_file("data/${fil}" "${CMAKE_CURRENT_BINARY_DIR}/data/${file}" COPYONLY)
        endforeach ()
        add_executable(test-${PROJECT_NAME}-HDF5Handle.exe testHDF5Handle.cpp)
        target_link_libraries(test-${PROJECT_NAME}-HDF5Handle.exe PUBLIC molpro::${PROJECT_NAME} gmock_main)
        target_compile_definitions(test-${PROJECT_NAME}-HDF5Handle.exe PUBLIC "ARRAY_DATA=\"${CMAKE_CURRENT_BINARY_DIR}/data\"")
        gtest_discover_tests(test-${PROJECT_NAME}-HDF5Handle.exe)
    endif ()
    add_executable(test-${PROJECT_NAME}-parallel
            main.cpp
            testutil.cpp
            testDistrArray.h
            )
    if (LINEARALGEBRA_ARRAY_MPI3)
        target_sources(test-${PROJECT_NAME}-parallel PRIVATE testDistrArrayMPI3.cpp)
    endif ()
    if (LINEARALGEBRA_ARRAY_GA)
        target_sources(test-${PROJECT_NAME}-parallel PRIVATE testDistrArrayGA.cpp)
    endif ()
    target_link_libraries(test-${PROJECT_NAME}-parallel PUBLIC molpro::${PROJECT_NAME} gmock_main)
    gtest_discover_tests(test-${PROJECT_NAME}-parallel)

    add_test(NAME test-${PROJECT_NAME}-parallel
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIOPTIONS} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            ${CMAKE_CURRENT_BINARY_DIR}/test-${PROJECT_NAME}-parallel ${MPIEXEC_POSTFLAGS}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            )
endif ()
