if (MPI)
    set(HDF5_PREFER_PARALLEL ON)
endif ()
if (LINEARALGEBRA_ARRAY_HDF5)
    find_package(HDF5 REQUIRED COMPONENTS C)
endif ()
find_package(Threads)

if (MPI AND LINEARALGEBRA_ARRAY_HDF5 AND NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR "Parallel HDF5 is required when building LinearAlgebra with MPI and HDF5 on")
endif ()

if (LINEARALGEBRA_ARRAY_HDF5)
    LibraryManager_Append(${PROJECT_NAME}
            SOURCES HDF5Handle.cpp
            PUBLIC_HEADER HDF5Handle.h
            )
endif ()

LibraryManager_Append(${PROJECT_NAME}
        PUBLIC_HEADER ArrayHandler.h ArrayHandlerIterable.h ArrayHandlerSparse.h
        ArrayHandlerIterableSparse.h ArrayHandlerDefault.h ArrayHandlerDDisk.h ArrayHandlerDDiskDistr.h
        ArrayHandlerDistrDDisk.h ArrayHandlerDDiskSparse.h Span.h type_traits.h
        default_handler.h
        )
if (MPI)
    LibraryManager_Append(${PROJECT_NAME} PUBLIC_HEADER ArrayHandlerDistr.h ArrayHandlerDistrSparse.h)
endif ()

if (LINEARALGEBRA_ARRAY_HDF5)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${HDF5_LIBRARIES})
    target_include_directories(${PROJECT_NAME} PUBLIC ${HDF5_INCLUDE_DIRS})
    target_compile_definitions(${PROJECT_NAME} PUBLIC ${HDF5_DEFINITIONS})
endif ()

    LibraryManager_Append(${PROJECT_NAME}
            SOURCES DistrArray.cpp DistrArrayDisk.cpp DistrArrayFile.cpp DistrArraySpan.cpp util.cpp
            PUBLIC_HEADER DistrArray.h DistrArrayDisk.h DistrArrayFile.h DistrArraySpan.h util.h
            )
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${PROJECT_NAME} PUBLIC "-pthread") # needed for DistrArrayDisk
        target_link_options(${PROJECT_NAME} PUBLIC "-pthread") # needed for DistrArrayDisk
    endif ()

if (MPI)
    if (LINEARALGEBRA_ARRAY_HDF5)
        LibraryManager_Append(${PROJECT_NAME}
                SOURCES PHDF5Handle.cpp
                PUBLIC_HEADER PHDF5Handle.h
                )
    endif ()

    if (LINEARALGEBRA_ARRAY_MPI3)
        message(STATUS "Building ${PROJECT_NAME} with MPI3 Array")
        target_compile_definitions(${PROJECT_NAME} PUBLIC LINEARALGEBRA_ARRAY_MPI3 LINEARALGEBRA_ARRAY_FILE)
        LibraryManager_Append(${PROJECT_NAME} SOURCES DistrArrayMPI3.cpp
                PUBLIC_HEADER DistrArrayMPI3.h ghc/filesystem.h)
    endif ()

    if (LINEARALGEBRA_ARRAY_GA)
        message(STATUS "Building ${PROJECT_NAME} with GA distributed array")
        target_compile_definitions(${PROJECT_NAME} PUBLIC LINEARALGEBRA_ARRAY_GA)
        LibraryManager_Append(${PROJECT_NAME} SOURCES DistrArrayGA.cpp PUBLIC_HEADER DistrArrayGA.h)
        target_link_libraries(${PROJECT_NAME} PUBLIC ga::ga)
    endif ()

    if (LINEARALGEBRA_ARRAY_HDF5)
        message(STATUS "Building ${PROJECT_NAME} with HDF5 Array")
        set(HDF5_PREFER_PARALLEL ON)
        find_package(HDF5 REQUIRED COMPONENTS C)
        target_compile_definitions(${PROJECT_NAME} PUBLIC LINEARALGEBRA_ARRAY_HDF5)
        target_link_libraries(${PROJECT_NAME} PUBLIC Threads::Threads)
        LibraryManager_Append(${PROJECT_NAME}
                SOURCES DistrArrayHDF5.cpp
                PUBLIC_HEADER DistrArrayHDF5.h)
    endif ()
endif ()

add_subdirectory(util)