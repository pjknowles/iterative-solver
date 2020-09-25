LibraryManager_Append(${PROJECT_NAME}
        SOURCES Logger.cpp
        PUBLIC_HEADER IterativeSolver.h IterativeSolverTemplate.h LinearEigensystemA.h
        Logger.h
        )
add_subdirectory(subspace)
