add_executable(solve-with-foreign-container solve-with-foreign-container.cpp)
target_link_libraries(solve-with-foreign-container PRIVATE molpro::LinearAlgebra)