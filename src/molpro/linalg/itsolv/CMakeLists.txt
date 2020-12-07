LibraryManager_Append(${PROJECT_NAME}
        SOURCES Logger.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystem.h Options.h LinearEquations.h
        Logger.h wrap.h wrap_util.h propose_rspace.h util.h DSpaceResetter.h LinearEigensystemOptions.h CastOptions.h
        LinearEquationsOptions.h
        )
add_subdirectory(subspace)
