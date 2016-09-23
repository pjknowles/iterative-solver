project(IterativeSolver)
cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)
enable_language(CXX)
if(CMAKE_COMPILER_IS_GNUCXX)
    add_compile_options(-Wall)
endif()
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    add_definitions(-DEIGEN_DONT_VECTORIZE)
endif()
aux_source_directory(. SRC_LIST)
set(SRC_LIST IterativeSolver.cpp IterativeSolverF.F90 SimpleParameterVector.cpp ParameterVector.cpp main.cpp Davidson.cpp Diis.cpp)
add_executable(${PROJECT_NAME}-test ${SRC_LIST})
#add_library (${PROJECT_NAME} ${SRC_LIST})
include_directories(Eigen )


find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(${PROJECT_NAME}-doc ALL
    DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
    )
file(GLOB headers *.h *.fh)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
DEPENDS ${SRC_LIST} ${headers} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)
