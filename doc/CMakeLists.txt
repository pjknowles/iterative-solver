# Documentation has to come last to make sure all dependencies have been added
include(DependencyManagerDocs)

foreach (f DoxygenLayout.xml)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${f} ${DEPENDENCYMANAGERDOCS_BASE_DIR}/${PROJECT_NAME}/${f} @ONLY)
endforeach ()

DependencyManagerDocs_Add(${PROJECT_NAME}
        FILES ../README.md
        DOC_URL "https://molpro.gitlab.io/IterativeSolver/${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
        DEPENDS Profiler # projects whose documentation has to be built before current project.
        )
