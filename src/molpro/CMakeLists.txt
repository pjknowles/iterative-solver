LibraryManager_Append(${PROJECT_NAME}
        SOURCES
        IterativeSolver.cpp
        PUBLIC_HEADER
        SimpleVector.h
        PagedVector.h
        OpaqueVector.h
        OutOfCoreArray.h
        IterativeSolver.h
        )
if (CMAKE_Fortran_COMPILER)
    LibraryManager_Append(${PROJECT_NAME}
            SOURCES IterativeSolverF.F90
            )
endif ()

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

include(DependencyManager)
DependencyManager_Populate(Profiler)
DependencyManager_Populate(utilities)
target_link_libraries(${PROJECT_NAME} PUBLIC molpro::Profiler molpro::utilities)

if (MPI_CXX_FOUND)
    message(STATUS "Building IterativeSolver on MPI")
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAVE_MPI_H)
    target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)
else ()
    message(STATUS "Building serial IterativeSolver")
endif ()
if (TARGET ppidd::ppidd)
    message(STATUS "Building IterativeSolver on PPIDD library")
    target_link_libraries(${PROJECT_NAME} PUBLIC ppidd::ppidd)
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAVE_PPIDD_H)
endif ()

LibraryManager_BLAS(${PROJECT_NAME} Intel10_64lp_seq Intel10_64lp)
message(STATUS "MKL_TYPE=${MKL_TYPE}")
if (MKL_TYPE MATCHES "^Intel10_64lp.*")
    message(STATUS "Eigen will use MKL")
    if(APPLE) # should be inside BLA_VENDOR TEST BELOW
        target_link_options(${PROJECT_NAME} PUBLIC "-Wl,-rpath,$ENV{MKLROOT}/lib")
    endif()
    target_compile_definitions(${PROJECT_NAME} PUBLIC EIGEN_USE_MKL_ALL)
endif ()
if (CMAKE_Fortran_COMPILER AND NOT LAPACK_FOUND)
    find_package(LAPACK REQUIRED)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES})
endif ()

if (NOT TARGET Eigen3::Eigen)
    message(VERBOSE "Eigen3::Eigen is not yet a target")
    find_package(Eigen3 3.3.0 NO_MODULE NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PACKAGE_REGISTRY QUIET)
    if (EIGEN3_FOUND)
        message(STATUS "Eigen3 found on system: ${EIGEN3_INCLUDE_DIR}")
    else ()
        if (NOT CMAKE_Fortran_COMPILER)
            message(FATAL_ERROR "Eigen3 cannot be configured from source with FORTRAN=OFF; install it in the system first")
        endif ()
        include(FetchContent)
        FetchContent_Declare(eigen3
                GIT_REPOSITORY https://gitlab.com/libeigen/eigen
                GIT_TAG 3.3.7
                )
        FetchContent_MakeAvailable(eigen3)
        message(STATUS "Eigen3::Eigen not found on system, and downloaded")
    endif ()
else ()
    message(STATUS "SymmetryMatrix: Eigen3::Eigen target already defined")
endif ()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    set(EIGEN_DONT_VECTORIZE ON)
endif ()

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_DONT_VECTORIZE)
endif ()

# Get all propreties that cmake supports
execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)

# Convert command output into a CMake list
STRING(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
STRING(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
# Fix https://stackoverflow.com/questions/32197663/how-can-i-remove-the-the-location-property-may-not-be-read-from-target-error-i
list(FILTER CMAKE_PROPERTY_LIST EXCLUDE REGEX "^LOCATION$|^LOCATION_|_LOCATION$")
# For some reason, "TYPE" shows up twice - others might too?
list(REMOVE_DUPLICATES CMAKE_PROPERTY_LIST)

# build whitelist by filtering down from CMAKE_PROPERTY_LIST in case cmake is
# a different version, and one of our hardcoded whitelisted properties
# doesn't exist!
unset(CMAKE_WHITELISTED_PROPERTY_LIST)
foreach(prop ${CMAKE_PROPERTY_LIST})
    if(prop MATCHES "^(INTERFACE|[_a-z]|IMPORTED_LIBNAME_|MAP_IMPORTED_CONFIG_)|^(COMPATIBLE_INTERFACE_(BOOL|NUMBER_MAX|NUMBER_MIN|STRING)|EXPORT_NAME|IMPORTED(_GLOBAL|_CONFIGURATIONS|_LIBNAME)?|NAME|TYPE|NO_SYSTEM_FROM_IMPORTED)$")
        list(APPEND CMAKE_WHITELISTED_PROPERTY_LIST ${prop})
    endif()
endforeach(prop)

function(print_properties)
    message ("CMAKE_PROPERTY_LIST = ${CMAKE_PROPERTY_LIST}")
endfunction(print_properties)

function(print_whitelisted_properties)
    message ("CMAKE_WHITELISTED_PROPERTY_LIST = ${CMAKE_WHITELISTED_PROPERTY_LIST}")
endfunction(print_whitelisted_properties)

function(print_target_properties tgt)
    if(NOT TARGET ${tgt})
        message("There is no target named '${tgt}'")
        return()
    endif()

    get_target_property(target_type ${tgt} TYPE)
    if(target_type STREQUAL "INTERFACE_LIBRARY")
        set(PROP_LIST ${CMAKE_WHITELISTED_PROPERTY_LIST})
    else()
        set(PROP_LIST ${CMAKE_PROPERTY_LIST})
    endif()

    foreach (prop ${PROP_LIST})
        string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" prop ${prop})
        # message ("Checking ${prop}")
        get_property(propval TARGET ${tgt} PROPERTY ${prop} SET)
        if (propval)
            get_target_property(propval ${tgt} ${prop})
            message ("${tgt} ${prop} = ${propval}")
        endif()
    endforeach(prop)
endfunction(print_target_properties)
print_target_properties(BLAS::BLAS)
#print_target_properties(Eigen3::Eigen)
print_target_properties(${PROJECT_NAME})