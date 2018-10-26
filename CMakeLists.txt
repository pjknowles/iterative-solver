project(IterativeSolver)
set(CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 2.8)

if ("${MPI_CXX_COMPILER}" STREQUAL "")
else ()
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    find_package(MPI)
    if (MPI_CXX_FOUND)
        add_definitions(-DHAVE_MPI_H)
    endif ()
endif ()

find_package(MPI )
if(MPI_CXX_FOUND)
include_directories(${MPI_CXX_INCLUDE_PATH})
add_definitions(-DHAVE_MPI_H)
# OpenMPI root guard: https://github.com/open-mpi/ompi/issues/4451
if("$ENV{USER}" STREQUAL "root" OR "$ENV{HOME}" STREQUAL "/root")
    set(MPIOPTIONS "--allow-run-as-root")
endif()
endif(MPI_CXX_FOUND)

include(ProcessorCount)
ProcessorCount(NPROC)
if(NPROC EQUAL 0)
    set(NPROC 1)
endif()
if(NPROC EQUAL 8)
    #    hack to avoid overload on workstation / CI
    set(NPROC 4)
endif()
message("Number of processors ${NPROC}")
set(CTEST_BUILD_FLAGS -j${NPROC})

enable_language(Fortran)
enable_language(CXX)
if (CMAKE_COMPILER_IS_GNUCXX)
    add_compile_options(-Wall)
endif ()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    add_definitions(-DEIGEN_DONT_VECTORIZE)
endif ()

aux_source_directory(IterativeSolver SRC_LIST)
set(SRC_LIST ${SRC_LIST} lib/LinearAlgebra.h lib/LAvector.h lib/PagedVector.h lib/IterativeSolverF.F90 lib/IterativeSolver.cpp lib/IterativeSolver.h)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Eigen3 3.3 REQUIRED)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})

add_library(${PROJECT_NAME} ${SRC_LIST})
target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/lib)

enable_testing()
add_executable(${PROJECT_NAME}-test ${PROJECT_NAME}-test.cpp lib/IterativeSolverF.F90)
target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
if(MPIEXEC)
    add_test(${PROJECT_NAME}-test ${MPIEXEC} ${MPIOPTIONS} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-test)
else(MPIEXEC)
    add_test(${PROJECT_NAME}-test ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-test ${fil})
endif(MPIEXEC)

add_executable(PagedVector-test PagedVector-test.cpp)
target_link_libraries(PagedVector-test ${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
if(MPIEXEC)
    add_test(PagedVector-test ${MPIEXEC} ${MPIOPTIONS} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${CMAKE_CURRENT_BINARY_DIR}/PagedVector-test)
else(MPIEXEC)
    add_test(PagedVector-test ${CMAKE_CURRENT_BINARY_DIR}/PagedVector-test)
endif(MPIEXEC)

file(GLOB EXAMPLES examples/*.cpp examples/*.F90)
foreach (fil ${EXAMPLES})
    get_filename_component(nam ${fil} NAME_WE)
    get_filename_component(ext ${fil} EXT)
    if (${ext} EQUAL ".F90")
        #        add_executable(${PROJECT_NAME}-${nam} IterativeSolverF.F90 ${fil})
        add_executable(${PROJECT_NAME}-${nam} ${fil})
    else ()
        add_executable(${PROJECT_NAME}-${nam} ${fil})
    endif ()
    set_target_properties(${PROJECT_NAME}-${nam} PROPERTIES COMPILE_FLAGS -UHAVE_MPI_H)
    target_link_libraries(${PROJECT_NAME}-${nam} ${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
    add_test(${nam} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-${nam})
endforeach ()


find_package(Doxygen)
if (DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(${PROJECT_NAME}-doc ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            )
    file(GLOB doxysrc *.h *.F90 examples/*)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            DEPENDS ${SRC_LIST} ${doxysrc} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif (DOXYGEN_FOUND)
