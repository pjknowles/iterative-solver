cmake_minimum_required(VERSION 3.11)

project(IterativeSolver LANGUAGES CXX VERSION 1.0)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(libconfig)

option(FORTRAN "Whether to build fortran binding" ON)
option(BENCHMARK "Whether to build benchmark job" OFF)
option(INTEGER8 "Whether to build for 64-bit fortran integers" ON)
if (Molpro_SOURCE_DIR)
    set(FORTRAN ON)
endif ()

if (BENCHMARK)
    set(FORTRAN ON)
endif ()

if (FORTRAN)
    enable_language(Fortran)
endif ()

########## Compiler and Build Flags ##########

if (NOT MPI_FOUND)
    # is the provided compiler an MPI wrapper, or is MPI specified?
    if (NOT MPI)
        set(MPI_SKIP_COMPILER_WRAPPER TRUE)
        set(MPI_SKIP_GUESSING TRUE)
    endif ()
    find_package(MPI)
endif (NOT MPI_FOUND)

if (MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set(HAVE_MPI_H ON)
    message("MPI build of ${CMAKE_CURRENT_SOURCE_DIR} for ${CMAKE_PROJECT_NAME} with CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} Fortran=${CMAKE_Fortran_COMPILER}")
else ()
    message("Serial build of ${CMAKE_CURRENT_SOURCE_DIR} for ${CMAKE_PROJECT_NAME} with CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} Fortran=${CMAKE_Fortran_COMPILER}")
endif (MPI_CXX_FOUND)


if (FORTRAN)
    message("Configuring Fortran support for ${CMAKE_CURRENT_SOURCE_DIR}")
    include(CheckFortranCompilerFlag)
    if (INTEGER8)
        foreach (f "-fdefault-integer-8" "-i8")
            CHECK_Fortran_COMPILER_FLAG(${f} _fortran_flags)
            if (_fortran_flags)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${f}")
            endif ()
            unset(_fortran_flags CACHE)
        endforeach ()
    endif ()
endif ()


########## Dependencies ##########

add_subdirectory(dependencies)

########## Targets ##########
add_subdirectory(lib)

enable_testing()
add_subdirectory(test)
add_subdirectory(example)


########## Documentation ##########

find_package(Doxygen)
if (DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(${PROJECT_NAME}-doc ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            )
    file(GLOB doxysrc lib/*.h lib/*.F90 example/*.cpp example/*.F90)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            DEPENDS ${SRC_LIST} ${doxysrc} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif (DOXYGEN_FOUND)
