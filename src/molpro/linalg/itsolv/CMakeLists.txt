LibraryManager_Append(${PROJECT_NAME}
        SOURCES util.cpp Logger.cpp Options.cpp LinearEigensystemOptions.cpp LinearEquationsOptions.cpp NonLinearEquationsDIISOptions.cpp OptimizeBFGSOptions.cpp OptimizeSDOptions.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystem.h Options.h LinearEquations.h
        SOURCES Logger.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystem.h OptimizeBFGS.h OptimizeSD.h Options.h LinearEquations.h
        Logger.h wrap.h wrap_util.h propose_rspace.h util.h DSpaceResetter.h LinearEigensystemOptions.h CastOptions.h
        LinearEquationsOptions.h SolverFactory.h SolverFactory-implementation.h options_map.h
        OptimizeBFGSOptions.h OptimizeSDOptions.h
        NonLinearEquationsDIISOptions.h
        )
add_subdirectory(subspace)
