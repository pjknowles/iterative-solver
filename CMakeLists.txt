cmake_minimum_required(VERSION 3.10)

project(IterativeSolver LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 11)


########## Compiler and Build Flags ##########

if ("${MPI_CXX_COMPILER}" STREQUAL "")
else ()
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    find_package(MPI)
    if (MPI_CXX_FOUND)
        add_definitions(-DHAVE_MPI_H)
    endif ()
endif ()

find_package(MPI)
if (MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    add_definitions(-DHAVE_MPI_H)
    # OpenMPI root guard: https://github.com/open-mpi/ompi/issues/4451
    if ("$ENV{USER}" STREQUAL "root" OR "$ENV{HOME}" STREQUAL "/root")
        set(MPIOPTIONS "--allow-run-as-root")
    endif ()
endif (MPI_CXX_FOUND)

if (CMAKE_COMPILER_IS_GNUCXX)
    add_compile_options(-Wall)
endif ()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    add_definitions(-DEIGEN_DONT_VECTORIZE)
endif ()

include(ProcessorCount)
ProcessorCount(NPROC)
if (NPROC EQUAL 0)
    set(NPROC 1)
endif ()
if (NPROC EQUAL 8)
    #    hack to avoid overload on workstation / CI
    set(NPROC 4)
endif ()
message("Number of processors ${NPROC}")
set(CTEST_BUILD_FLAGS -j${NPROC})


########## Dependencies ##########

if (NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.3)
    if (NOT TARGET Eigen3::Eigen)
        include(FetchContent)
        FetchContent_Declare(
                eigen3_dep
                SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/Eigen3"
                GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
                GIT_TAG "3.3.0"
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
    endif ()
endif ()


########## Targets ##########
add_subdirectory(lib)

enable_testing()
add_subdirectory(test)
add_subdirectory(examples)


########## Documentation ##########

find_package(Doxygen)
if (DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(${PROJECT_NAME}-doc ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            )
    file(GLOB doxysrc *.h *.F90 examples/*)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            DEPENDS ${SRC_LIST} ${doxysrc} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif (DOXYGEN_FOUND)
