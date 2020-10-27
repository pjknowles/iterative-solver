LibraryManager_Append(${PROJECT_NAME}
        SOURCES Logger.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystemA.h
        Logger.h wrap.h propose_rspace.h reset_dspace.h util.h
        )
add_subdirectory(subspace)
