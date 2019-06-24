set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)
add_library(${PROJECT_NAME} STATIC
        SimpleVector.h
        PagedVector.h
        OpaqueVector.h
        IterativeSolver.cpp IterativeSolver.h
        ${CMAKE_CURRENT_BINARY_DIR}/IterativeSolver-config.h
        )
if (FORTRAN)
    target_sources(${PROJECT_NAME} PRIVATE IterativeSolverF.F90)
endif ()
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)
if (Molpro_SOURCE_DIR)
    if (NOT TARGET ppidd::ppidd)
        if (NOT ppidd_dep_POPULATED)
            FetchContent_Populate(ppidd_dep)
            add_subdirectory(${ppidd_dep_SOURCE_DIR} ${ppidd_dep_BINARY_DIR} EXCLUDE_FROM_ALL)
        endif ()
    endif ()
    list(APPEND DEPENDENCIES ppidd)
    set(DEPENDENCY_ppidd 1.0)
    target_link_libraries(${PROJECT_NAME} PRIVATE ppidd::ppidd)
endif ()
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "SimpleVector.h;PagedVector.h;OpaqueVector.h;IterativeSolver.h;${CMAKE_CURRENT_BINARY_DIR}/IterativeSolver-config.h")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
set (EIGEN_DONT_VECTORIZE ON)
endif ()

include (libconfig.cmake)
