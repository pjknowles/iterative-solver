LibraryManager_Append(${PROJECT_NAME}
        SOURCES ArrayGA.cpp
        PUBLIC_HEADER Array.h ArrayGA.h
        )

if (GCI_ARRAY_MPI3)
    message(STATUS "Building gci with MPI3 Array")
    target_compile_definitions(${PROJECT_NAME} PUBLIC GCI_ARRAY_ArrayMPI3)
elseif (GCI_ARRAY_GA)
    message(STATUS "Building gci with GA SharedCounter")
    target_compile_definitions(${PROJECT_NAME} PUBLIC GCI_ARRAY_ArrayGA)
else ()
    message(FATAL_ERROR "Array type not specified, turn on one of the options: GCI_ARRAY_MPI3, GCI_ARRAY_GA")
endif ()
