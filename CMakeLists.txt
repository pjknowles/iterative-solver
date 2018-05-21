project(IterativeSolver)
set(CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 2.8)

if("${MPI_CXX_COMPILER}" STREQUAL "")
else()
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    find_package(MPI)
    if(MPI_CXX_FOUND)
        add_definitions(-DUSE_MPI)
    endif()
endif()

enable_language(Fortran)
enable_language(CXX)
if(CMAKE_COMPILER_IS_GNUCXX)
    add_compile_options(-Wall)
endif()
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    add_definitions(-DEIGEN_DONT_VECTORIZE)
endif()

aux_source_directory(IterativeSolver SRC_LIST)
set(SRC_LIST ${SRC_LIST} LinearAlgebra.h LAvector.h IterativeSolverF.F90 IterativeSolver.cpp)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package (Eigen3 3.3 REQUIRED)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})

add_library (${PROJECT_NAME} ${SRC_LIST})

enable_testing()
add_executable(${PROJECT_NAME}-test ${PROJECT_NAME}-test.cpp IterativeSolverF.F90)
target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
add_test(${PROJECT_NAME}-test ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-test ${fil})

add_executable(PagedVector-test PagedVector-test.cpp)
target_link_libraries(PagedVector-test ${PROJECT_NAME})
add_test(PagedVector-test ${CMAKE_CURRENT_BINARY_DIR}/PagedVector-test)

file (GLOB TESTS examples/*.cpp)
foreach(fil ${TESTS})
  get_filename_component(nam ${fil} NAME_WE)
  add_executable(${PROJECT_NAME}-${nam} examples/${nam}.cpp)
  target_link_libraries(${PROJECT_NAME}-${nam} ${PROJECT_NAME})
  add_test(${nam} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-${nam})
endforeach()


find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(${PROJECT_NAME}-doc ALL
    DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
    )
file(GLOB doxysrc IterativeSolver/*.h IterativeSolver/*.fh examples/*)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
DEPENDS ${SRC_LIST} ${doxysrc} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)
