project(diiscxx)
cmake_minimum_required(VERSION 2.8)
add_subdirectory(Eigen)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
include_directories(Eigen )

