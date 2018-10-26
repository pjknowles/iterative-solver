########## C++ Library ##########
add_library(${PROJECT_NAME} STATIC
        LinearAlgebra.h
        LAvector.h
        PagedVector.h
        IterativeSolver.cpp IterativeSolver.h
        )
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

########## Fortran Library ##########
add_library(${PROJECT_NAME}F STATIC IterativeSolverF.F90)
target_link_libraries(${PROJECT_NAME}F PUBLIC ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
