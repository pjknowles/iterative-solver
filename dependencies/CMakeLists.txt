#declare_dependency(googletest https://github.com/google/googletest.git)
#declare_dependency(eigen3 https://github.com/eigenteam/eigen-git-mirror.git)

include(FetchContent)
FetchContent_Declare(
        dependency_manager
        GIT_REPOSITORY https://gitlab.com/dependencymanager/dependency-manager.git
        GIT_TAG 0.1.5
)
FetchContent_MakeAvailable(dependency_manager)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

DependencyManager_Declare(ppidd https://gitlab.com/molpro/ppidd.git)

DependencyManager_Declare(Profiler https://gitlab.com/molpro/profiler.git)

DependencyManager_Declare( utilities https://gitlab+deploy-token-181030:whWoM4qkenzjGjWsFTxj@gitlab.com/molpro/utilities.git )
