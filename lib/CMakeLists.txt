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

set(DEPENDENCY_Eigen3::Eigen 3.3.0)
if (NOT TARGET Eigen3::Eigen)
    get_dependency_name("eigen3")
    FetchContent_GetProperties(${_dependency_name})
    if (NOT ${_dependency_name}_POPULATED)
        find_package(Eigen3 3.3.0 NO_MODULE)
    endif ()
    if (NOT TARGET Eigen3::Eigen)
        get_dependency(eigen3 NO)
        add_library(eigen INTERFACE)
        target_include_directories(eigen INTERFACE
                $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/dependencies/Eigen3>
                $<INSTALL_INTERFACE:include>
                )
        add_library(Eigen3::Eigen ALIAS eigen)
        list(APPEND DEPENDENCIES eigen)
        message(STATUS "IterativeSolver: Eigen3 library obtained externally")
    endif ()
else ()
    message(STATUS "IterativeSolver: Eigen3::Eigen target already defined")
endif ()

if (Molpro_SOURCE_DIR)
    get_dependency_name(ppidd)
    list(APPEND DEPENDENCIES ppidd)
    set(DEPENDENCY_ppidd 1.0)
    target_link_libraries(${PROJECT_NAME} PRIVATE ppidd::ppidd)
endif ()
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "SimpleVector.h;PagedVector.h;OpaqueVector.h;IterativeSolver.h;${CMAKE_CURRENT_BINARY_DIR}/IterativeSolver-config.h")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
set (EIGEN_DONT_VECTORIZE ON)
endif ()

configure_library(${PROJECT_NAME} "${DEPENDENCIES}")
