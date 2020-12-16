# Documentation has to come last to make sure all dependencies have been added
include(DependencyManagerDocs)

foreach (f DoxygenLayout.xml)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${f} ${DEPENDENCYMANAGERDOCS_BASE_DIR}/${PROJECT_NAME}/${f} @ONLY)
endforeach ()

foreach (f images/solvers.png images/dependency_inversion.png images/options.png images/solver_factory.png)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${f} ${DEPENDENCYMANAGERDOCS_BASE_DIR}/${PROJECT_NAME}/html/${f} COPYONLY)
endforeach ()

DependencyManagerDocs_Add(${PROJECT_NAME}
        FILES ../README.md LibraryDesign.md IterativeSolver.dox Arrays.dox ArrayHandlers.dox VecRef.dox
        DOC_URL "https://molpro.gitlab.io/IterativeSolver/${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
        DEPENDS Profiler # projects whose documentation has to be built before current project.
        )
