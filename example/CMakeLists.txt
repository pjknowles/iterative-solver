if (CMAKE_Fortran_COMPILER_ID
    AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Intel"  # avoid the challenge of linking C++ and Fortran with Intel
)
    set(ITERATIVESOLVER_FORTRAN_EXAMPLES
            OptimizeExampleF.F90
            DIISExampleF.F90
            LinearEigensystemExampleF.F90
            LinearEigensystemExampleF-Pspace.F90
            LinearEigensystemExampleF-Pspace-adaptive.F90
            LinearEquationsExampleF.F90
            )
endif ()
set(EXAMPLES DIISexample.cpp
        OptimizeExample.cpp
        LinearEigensystemExample.cpp LinearEigensystemExample-paged.cpp
        LinearEquationsExample.cpp
        ${ITERATIVESOLVER_FORTRAN_EXAMPLES}
        )

foreach (fil ${EXAMPLES})
    get_filename_component(nam ${fil} NAME_WE)
    set(test_name ${PROJECT_NAME}-${nam})
    add_executable(${test_name} ${fil} dummy.cpp)
    set_target_properties(${test_name}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
            )
    target_link_libraries(${test_name} PUBLIC ${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        target_link_libraries(${test_name} PUBLIC ifcore)
    endif()
    add_test(${nam} ${PROJECT_BINARY_DIR}/${test_name})
endforeach ()

if (BENCHMARK)
    add_executable(LinearEigensystemBenchmarkF LinearEigensystemBenchmarkF.F90)
    set_target_properties(LinearEigensystemBenchmarkF
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
            )
    get_dependency(profiler)
    if (MPI_CXX_FOUND)
        target_link_libraries(LinearEigensystemBenchmarkF PUBLIC ${PROJECT_NAME} Profiler::mpi)
    else ()
        target_link_libraries(LinearEigensystemBenchmarkF PUBLIC ${PROJECT_NAME} Profiler::serial)
    endif ()
endif ()
