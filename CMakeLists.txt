cmake_minimum_required(VERSION 3.11)

project(IterativeSolver LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 14)


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
    add_definitions(-DHAVE_MPI_H)
    message("MPI build of ${PROJECT_NAME}")
else ()
    message("Serial build of ${PROJECT_NAME}")
endif (MPI_CXX_FOUND)

if (CMAKE_COMPILER_IS_GNUCXX)
    add_compile_options(-Wall)
endif ()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    add_definitions(-DEIGEN_DONT_VECTORIZE)
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

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    enable_testing()
    add_subdirectory(test)
    add_subdirectory(example)
endif ()


########## Documentation ##########

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
        add_custom_target(${PROJECT_NAME}-doc ALL
                DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
                )
        file(GLOB doxysrc *.h *.F90 example/*)
        add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
                COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                DEPENDS ${SRC_LIST} ${doxysrc} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating API documentation with Doxygen" VERBATIM
                )
    endif (DOXYGEN_FOUND)
endif ()
