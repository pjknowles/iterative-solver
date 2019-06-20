cmake_minimum_required(VERSION 3.11)

project(IterativeSolver LANGUAGES CXX VERSION 1.0)
string(TOUPPER ${PROJECT_NAME} PROJECT_UPPER_NAME)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(${PROJECT_UPPER_NAME}_FORTRAN "Whether to build fortran binding" ON)
    option(${PROJECT_UPPER_NAME}_I8 "Whether to build for 64-bit fortran integers" ON)
endif ()
if (${PROJECT_UPPER_NAME}_FORTRAN)
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

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-Wall" _Wall)
if (_Wall)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
endif ()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    set(EIGEN_DONT_VECTORIZE ON)
endif ()

if (CMAKE_Fortran_COMPILER_ID)
    message("Configuring Fortran support for ${CMAKE_CURRENT_SOURCE_DIR}")
    include(CheckFortranCompilerFlag)
    if (${PROJECT_UPPER_NAME}_I8)
        foreach (f "-fdefault-integer-8" "-i8")
            CHECK_Fortran_COMPILER_FLAG(${f} _fortran_flags)
            if (_fortran_flags)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${f}")
            endif ()
            unset(_fortran_flags CACHE)
        endforeach ()
    endif ()
    CHECK_Fortran_COMPILER_FLAG("-Wall" _Wallf)
    if (_Wallf)
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall")
    endif ()
endif ()


########## Dependencies ##########

if (NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.3.5)
    if (NOT TARGET Eigen3::Eigen)
        include(FetchContent)
        FetchContent_Declare(
                eigen3_dep
                SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/Eigen3"
                GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
                GIT_TAG "3.3.5"
                GIT_SHALLOW TRUE
        )
        FetchContent_GetProperties(eigen3_dep)
        if (NOT eigen3_dep_POPULATED)
            FetchContent_Populate(eigen3_dep)
            add_library(ExternalEigen3 INTERFACE)
            target_include_directories(ExternalEigen3 INTERFACE ${CMAKE_SOURCE_DIR}/dependencies/Eigen3)
            add_library(Eigen3::Eigen ALIAS ExternalEigen3)
            message("Using external Eigen3 library")
        endif ()
    else ()
        message("Using existing Eigen3 library")
    endif ()
endif ()


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
