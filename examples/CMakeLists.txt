if (FORTRAN
        AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Intel"  # avoid the challenge of linking C++ and Fortran with Intel
        AND OFF # TODO reinstate Fortran examples when Fortran support is reinstated
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
        LinearEigensystemExample.cpp LinearEigensystemExample-serial.cpp LinearEigensystemExample-Pspace-fixed-serial.cpp
        LinearEquationsExample.cpp LinearEquationsExample-Pspace.cpp
        ${ITERATIVESOLVER_FORTRAN_EXAMPLES}
        )
if (LINEARALGEBRA_ARRAY_MPI3 AND MPI)
    list(APPEND EXAMPLES
            LinearEigensystemExample-mpi.cpp
            )
    if (FORTRAN)
        list(APPEND EXAMPLES
                LinearEigensystemExampleF-mpi.F90
                LinearEigensystemExampleF-Pspace-mpi.F90
                )
    endif ()
endif ()

if (FORTRAN AND MPI)
    find_package(MPI COMPONENTS Fortran)
    message("MPI_Fortran_LIBRARIES=${MPI_Fortran_LIBRARIES}")
endif ()
foreach (fil ${EXAMPLES})
    get_filename_component(nam ${fil} NAME_WE)
    set(test_name ${PROJECT_NAME}-${nam})
    add_executable(${test_name} ${fil} dummy.cpp)
    set_target_properties(${test_name}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
            )
    target_link_libraries(${test_name} PUBLIC molpro::${PROJECT_NAME} ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})
    if (FORTRAN AND CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        target_link_libraries(${test_name} PUBLIC ifcore)
    endif ()


    if (MPI)
        # in principle this should be configured at runtime but it seems cmake/ctest support is not yet there to do it simply
        if (NOT DEFINED ENV{MPIPROCS})
            set(ENV{MPIPROCS} 1)
        endif ()
        message(DEBUG "add_test(${nam} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} $ENV{MPIPROCS} ${MPIEXEC_PREFLAGS} ${PROJECT_BINARY_DIR}/${test_name} ${MPIEXEC_POSTFLAGS})")
        add_test(${nam} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} $ENV{MPIPROCS} ${MPIEXEC_PREFLAGS} ${PROJECT_BINARY_DIR}/${test_name} ${MPIEXEC_POSTFLAGS})
    else ()
        add_test(${nam} ${PROJECT_BINARY_DIR}/${test_name})
    endif ()
endforeach ()

if (BENCHMARK AND FORTRAN)
    add_executable(LinearEigensystemBenchmarkF LinearEigensystemBenchmarkF.F90)
    set_target_properties(LinearEigensystemBenchmarkF
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
            )
    if (MPI_CXX_FOUND)
        target_link_libraries(LinearEigensystemBenchmarkF PUBLIC molpro::${PROJECT_NAME})# Profiler::mpi)
    else ()
        target_link_libraries(LinearEigensystemBenchmarkF PUBLIC molpro::${PROJECT_NAME})# Profiler::serial)
    endif ()
endif ()

file(COPY bh.hamiltonian DESTINATION .)
file(COPY hf.hamiltonian DESTINATION .)
