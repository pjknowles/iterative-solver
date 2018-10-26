
set(EXAMPLES DIISexample.cpp DIISExampleF.F90
        LinearEigensystemExample.cpp LinearEigensystemExample-paged.cpp LinearEigensystemExampleF.F90
        LinearEigensystemExampleF-Pspace.F90 LinearEigensystemExampleF-Pspace-adaptive.F90
        LinearEquationsExample.cpp LinearEquationsExampleF.F90
        )

foreach (fil ${EXAMPLES})
    get_filename_component(nam ${fil} NAME_WE)
    set(test_name ${PROJECT_NAME}-${nam})
    add_executable(${test_name} ${fil})
    set_target_properties(${test_name}
            PROPERTIES
            COMPILE_FLAGS -UHAVE_MPI_H
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
            )
    get_filename_component(ext ${fil} EXT)
    if (${ext} STREQUAL ".F90")
        target_link_libraries(${test_name} PUBLIC ${PROJECT_NAME}F ${MPI_CXX_LIBRARIES})
    else ()
        target_link_libraries(${test_name} PUBLIC ${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
    endif ()
    add_test(${nam} ${PROJECT_BINARY_DIR}/${test_name})
endforeach ()

