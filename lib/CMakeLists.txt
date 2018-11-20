if (CMAKE_Fortran_COMPILER_ID)
    set(ITERATIVESOLVER_FORTRAN_FILES "IterativeSolverF.F90")
endif ()
add_library(${PROJECT_NAME} STATIC
        SimpleVector.h
        PagedVector.h
        OpaqueVector.h
        IterativeSolver.cpp IterativeSolver.h
        ${ITERATIVESOLVER_FORTRAN_FILES}
        )
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
### For including C++ headers
target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
### For linking Fortran library
target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
