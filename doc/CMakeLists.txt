# Documentation has to come last to make sure all dependencies have been added
include(DependencyManagerDocs)
DependencyManagerDocs_Add(${PROJECT_NAME}
        FILES README.md
        DOC_URL "https://molpro.gitlab.io/IterativeSolver/${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
        DEPENDS Profiler # projects whose documentation has to be built before current project.
        )
