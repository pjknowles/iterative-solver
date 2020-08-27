if (MPI)
    set(HDF5_PREFER_PARALLEL ON)
endif ()
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Threads)

if (MPI AND NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR "Parallel HDF5 is required when building LinearAlgebra with MPI on")
endif ()

LibraryManager_Append(${PROJECT_NAME}
        SOURCES HDF5Handle.cpp
        PUBLIC_HEADER HDF5Handle.h ArrayHandler.h ArrayHandlerIterable.h ArrayHandlerSparse.h
        ArrayHandlerIterableSparse.h ArrayHandlerDefault.h ArrayHandlerDDisk.h Span.h type_traits.h
        )
if (MPI)
    LibraryManager_Append(${PROJECT_NAME} PUBLIC_HEADER ArrayHandlerDistr.h ArrayHandlerDistrSparse.h)
endif ()
target_link_libraries(${PROJECT_NAME} PUBLIC ${HDF5_LIBRARIES})
target_include_directories(${PROJECT_NAME} PUBLIC ${HDF5_INCLUDE_DIRS})
target_compile_definitions(${PROJECT_NAME} PUBLIC ${HDF5_DEFINITIONS})

if (MPI)
    LibraryManager_Append(${PROJECT_NAME}
            SOURCES DistrArray.cpp DistrArrayDisk.cpp util.cpp PHDF5Handle.cpp
            PUBLIC_HEADER DistrArray.h DistrArrayDisk.h util.h PHDF5Handle.h
            )

    if (LINEARALGEBRA_ARRAY_MPI3)
        message(STATUS "Building ${PROJECT_NAME} with MPI3 Array")
        target_compile_definitions(${PROJECT_NAME} PUBLIC LINEARALGEBRA_ARRAY_MPI3)
        LibraryManager_Append(${PROJECT_NAME} SOURCES DistrArrayMPI3.cpp PUBLIC_HEADER DistrArrayMPI3.h)
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