LibraryManager_Append(${PROJECT_NAME}
        SOURCES Logger.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystem.h
        Logger.h wrap.h wrap_util.h propose_rspace.h util.h DSpaceResetter.h
        )
add_subdirectory(subspace)
