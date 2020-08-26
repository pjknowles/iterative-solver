LibraryManager_Append(${PROJECT_NAME}
        PUBLIC_HEADER Distribution.h select_max_dot.h)
if (MPI)
    LibraryManager_Append(${PROJECT_NAME}
            SOURCES DistrFlags.cpp
            PUBLIC_HEADER DistrFlags.h gather_all.h)
endif ()
