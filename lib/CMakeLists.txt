unset(DEPENDENCIES)
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
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)

if (NOT TARGET Eigen3::Eigen)
    message(STATUS "Eigen3::Eigen is not yet a target")
    find_package(Eigen3 3.3.0 NO_MODULE NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PACKAGE_REGISTRY QUIET)
    if (EIGEN3_FOUND)
        message(STATUS "Eigen3 found on system: ${EIGEN3_INCLUDE_DIR}")
        add_library(eigen INTERFACE)
        target_include_directories(eigen INTERFACE
                $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
                $<INSTALL_INTERFACE:${EIGEN3_INCLUDE_DIR}>
                )
        add_library(Eigen3::Eigen ALIAS eigen)
    else ()
        message(STATUS "Eigen3::Eigen not found on system, and downloaded")
        get_dependency(eigen3 OFF)
        get_dependency_name("eigen3")
        FetchContent_GetProperties(${_dependency_name})
        add_library(eigen INTERFACE)
        set(EIGEN3_INCLUDE_DIR ${${_dependency_name}_SOURCE_DIR})
        target_include_directories(eigen INTERFACE
                $<BUILD_INTERFACE:${${_dependency_name}_SOURCE_DIR}>
                $<INSTALL_INTERFACE:include>
                )
    endif ()
    add_library(Eigen3::Eigen ALIAS eigen)
else ()
    message(STATUS "IterativeSolver: Eigen3::Eigen target already defined")
endif ()
get_target_property(incl eigen INTERFACE_INCLUDE_DIRECTORIES)
install(TARGETS eigen EXPORT ${PROJECT_NAME}Targets LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        PUBLIC_HEADER DESTINATION include
        )
set(DEPENDENCY_Eigen3::Eigen 3.3.0)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

if (Molpro_SOURCE_DIR)
    message(STATUS "IterativeSolver Molpro configuration")
    get_dependency(ppidd)
    if (TARGET ppidd::ppidd)
        message(STATUS "ppidd::ppidd is a target")
    else ()
        message(STATUS "ppidd::ppidd is not a target")
endif()
if (TARGET ppidd)
message(STATUS "ppidd is a target")
else()
message(STATUS "ppidd is not a target")
endif()
    list(APPEND DEPENDENCIES ppidd)
    set(DEPENDENCY_ppidd 1.0)
    target_link_libraries(${PROJECT_NAME} PUBLIC ppidd::ppidd)
endif ()
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "SimpleVector.h;PagedVector.h;OpaqueVector.h;IterativeSolver.h;${CMAKE_CURRENT_BINARY_DIR}/IterativeSolver-config.h")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
set (EIGEN_DONT_VECTORIZE ON)
endif ()

configure_library(${PROJECT_NAME} "${DEPENDENCIES}")
