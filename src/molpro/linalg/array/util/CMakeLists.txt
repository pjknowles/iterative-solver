LibraryManager_Append(${PROJECT_NAME}
        SOURCES temp_file.cpp
        PUBLIC_HEADER Distribution.h select_max_dot.h temp_file.h temp_hdf5_handle.h TempHandle.h)
if (MPI)
    LibraryManager_Append(${PROJECT_NAME}
            SOURCES DistrFlags.cpp
            PUBLIC_HEADER DistrFlags.h)
    if (LINEARALGEBRA_ARRAY_HDF5)
        LibraryManager_Append(${PROJECT_NAME}
                PUBLIC_HEADER temp_phdf5_handle.h)
    endif ()
endif ()
