LibraryManager_Append(${PROJECT_NAME}
        SOURCES DistrArray.cpp util.cpp
        PUBLIC_HEADER DistrArray.h util.h
        )

if (LINEARALGEBRA_ARRAY_MPI3)
    message(STATUS "Building ${PROJECT_NAME} with MPI3 Array")
    target_compile_definitions(${PROJECT_NAME} PUBLIC LINEARALGEBRA_ARRAY_ArrayMPI3)
    LibraryManager_Append(${PROJECT_NAME} SOURCES DistrArrayMPI3.cpp PUBLIC_HEADER DistrArrayMPI3.h)
endif ()

if (LINEARALGEBRA_ARRAY_GA)
    message(STATUS "Building ${PROJECT_NAME} with GA distributed array")
    target_compile_definitions(${PROJECT_NAME} PUBLIC LINEARALGEBRA_ARRAY_ArrayGA)
    LibraryManager_Append(${PROJECT_NAME} SOURCES DistrArrayGA.cpp PUBLIC_HEADER DistrArrayGA.h)
    target_link_libraries(${PROJECT_NAME} PUBLIC ga::ga)
endif ()
