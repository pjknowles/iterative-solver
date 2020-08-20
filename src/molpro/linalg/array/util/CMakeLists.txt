LibraryManager_Append(${PROJECT_NAME}
        PUBLIC_HEADER Distribution.h)
if (MPI)
    LibraryManager_Append(${PROJECT_NAME}
            SOURCES DistrFlags.cpp
            PUBLIC_HEADER DistrFlags.h)
endif ()
