LibraryManager_Add(${PROJECT_NAME} NAMESPACE molpro)

# MPI has to be called in the same scope as LibraryManager_Install, if pkg-config file is needed
# Otherwise MPI::MPI_CXX target is not available.
if (MPI)
    find_package(MPI COMPONENTS CXX)
endif ()

if (PPIDD)
    if (CMAKE_Fortran_COMPILER)
        set(PPIDD_FORTRAN ON)
    else ()
        set(PPIDD_FORTRAN OFF)
    endif ()
    DependencyManager_Populate(ppidd)
endif ()

# Its source files are conventionally in directories named by C++ namespace
add_subdirectory(molpro)

#TODO reinstate these
#LibraryManager_Install(${PROJECT_NAME} PKG_CONFIG)
#LibraryManager_Export(${PROJECT_NAME})
