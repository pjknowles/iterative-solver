project(diiscxx)
cmake_minimum_required(VERSION 2.8)
add_subdirectory(Eigen)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
include_directories(Eigen )
add_definitions(-DEIGEN_USE_MKL_ALL)

set(MKL_ROOT /opt/intel/mkl)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(MKL)
if (MKL_FOUND)
include_directories(${MKL_INCLUDE_DIRS})
message("MKL link libraries are ${MKL_LIBRARIES}")
target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES})
else (MKL_FOUND)
find_package(BLAS REQUIRED)
if (BLAS_FOUND)
include_directories(${BLAS_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${BLAS_LIBRARY} -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl)
message("Blas link library is ${BLAS_LIBRARIES}")
else (BLAS_FOUND)
target_link_libraries(${PROJECT_NAME} -llapack -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl)
endif (BLAS_FOUND)
endif (MKL_FOUND)

include_directories(/opt/intel/mkl/include)

find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(${PROJECT_NAME}-doc ALL
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)
