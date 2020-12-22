LibraryManager_Append(${PROJECT_NAME}
        SOURCES util.cpp Logger.cpp Options.cpp LinearEigensystemOptions.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystem.h Options.h LinearEquations.h
        Logger.h wrap.h wrap_util.h propose_rspace.h util.h DSpaceResetter.h LinearEigensystemOptions.h CastOptions.h
        LinearEquationsOptions.h SolverFactory.h SolverFactory-implementation.h
        )
add_subdirectory(subspace)
