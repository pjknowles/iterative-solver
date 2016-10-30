project(IterativeSolver)
cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)
enable_language(CXX)
if(CMAKE_COMPILER_IS_GNUCXX)
    add_compile_options(-Wall)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    add_definitions(-DEIGEN_DONT_VECTORIZE)
endif()
aux_source_directory(IterativeSolver SRC_LIST)
set(SRC_LIST ${SRC_LIST} IterativeSolver/ParameterVector.h IterativeSolver/IterativeSolverF.F90)
include_directories(Eigen .)
add_library (${PROJECT_NAME} ${SRC_LIST})

add_executable(${PROJECT_NAME}-test main.cpp)
target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})

add_executable(${PROJECT_NAME}-DIISexample examples/DIISexample.cpp)
target_link_libraries(${PROJECT_NAME}-DIISexample ${PROJECT_NAME})

add_executable(${PROJECT_NAME}-RSPTexample examples/RSPTexample.cpp)
target_link_libraries(${PROJECT_NAME}-RSPTexample ${PROJECT_NAME})

add_executable(${PROJECT_NAME}-DavidsonExample examples/DavidsonExample.cpp)
target_link_libraries(${PROJECT_NAME}-DavidsonExample ${PROJECT_NAME})


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
