LibraryManager_Append(${PROJECT_NAME}
        SOURCES util.cpp Logger.cpp Options.cpp LinearEigensystemOptions.cpp LinearEquationsOptions.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystem.h Options.h LinearEquations.h
        SOURCES Logger.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystem.h Optimize.h Options.h LinearEquations.h
        Logger.h wrap.h wrap_util.h propose_rspace.h util.h DSpaceResetter.h LinearEigensystemOptions.h CastOptions.h
        LinearEquationsOptions.h SolverFactory.h SolverFactory-implementation.h options_map.h
        OptimizeOptions.h
        )
add_subdirectory(subspace)
