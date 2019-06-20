configure_file(IterativeSolver-config.h.in IterativeSolver-config.h)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)
add_library(${PROJECT_NAME} STATIC
        SimpleVector.h
        PagedVector.h
        OpaqueVector.h
        IterativeSolver.cpp IterativeSolver.h
#        IterativeSolverF.F90
        ${CMAKE_CURRENT_BINARY_DIR}/IterativeSolver-config.h
        )
if (CMAKE_Fortran_COMPILER_ID)
    target_sources(${PROJECT_NAME} PRIVATE IterativeSolverF.F90)
endif ()
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "SimpleVector.h;PagedVector.h;OpaqueVector.h;IterativeSolver.h;${CMAKE_CURRENT_BINARY_DIR}/IterativeSolver-config.h")

include(libconfig.cmake)
