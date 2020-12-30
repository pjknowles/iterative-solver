LibraryManager_Append(${PROJECT_NAME}
        SOURCES util.cpp Logger.cpp Options.cpp LinearEigensystemDavidsonOptions.cpp LinearEquationsDavidsonOptions.cpp NonLinearEquationsDIISOptions.cpp OptimizeBFGSOptions.cpp OptimizeSDOptions.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystemDavidson.h Options.h LinearEquationsDavidson.h
        SOURCES Logger.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystemDavidson.h OptimizeBFGS.h OptimizeSD.h Options.h LinearEquationsDavidson.h
        Logger.h wrap.h wrap_util.h propose_rspace.h util.h DSpaceResetter.h LinearEigensystemDavidsonOptions.h CastOptions.h
        LinearEquationsDavidsonOptions.h SolverFactory.h SolverFactory-implementation.h options_map.h
        OptimizeBFGSOptions.h OptimizeSDOptions.h
        NonLinearEquationsDIISOptions.h
        )
add_subdirectory(subspace)
