#declare_dependency(googletest https://github.com/google/googletest.git)
#declare_dependency(eigen3 https://github.com/eigenteam/eigen-git-mirror.git)

include(FetchContent)
FetchContent_Declare(
        dependency_manager
        GIT_REPOSITORY https://gitlab.com/Molpro/dependency-manager.git
        GIT_TAG 0.2.0
)
FetchContent_MakeAvailable(dependency_manager)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

DependencyManager_Declare(ppidd https://gitlab.com/molpro/ppidd.git)

DependencyManager_Declare(Profiler https://gitlab.com/molpro/profiler.git)

DependencyManager_Declare(utilities https://gitlab.com/molpro/utilities.git)

if (NOT PPIDD AND (MPI AND LINEARALGEBRA_ARRAY_GA))
    set(GA_RUNTIME "MPI_RMA" CACHE STRING "overwriting default for GA runtime protocol")
    if (NOT TARGET ga::ga)
        find_program(ga-config ga-config)
        if (ga-config)
            message("Local GA found")
            foreach (f cflags cppflags ldflags libs network_ldflags network_libs)
                execute_process(COMMAND ${ga-config} --${f} OUTPUT_VARIABLE ga-${f})
                string(STRIP "${ga-${f}}" ga-${f})
                separate_arguments(ga-${f} NATIVE_COMMAND "${ga-${f}}")
            endforeach ()
            string(FIND "${ga-libs}" "-framework;Accelerate" ga-accelerate)
            if (NOT ga-accelerate EQUAL -1)
                string(REPLACE ";-framework;Accelerate" "" ga-libs "${ga-libs}")
                string(APPEND ga-ldflags ";-framework;Accelerate")
            endif ()
            add_library(ga::ga INTERFACE IMPORTED GLOBAL)
            target_link_libraries(ga::ga INTERFACE "${ga-libs};${ga-network_libs}")
            target_link_options(ga::ga INTERFACE "${ga-ldflags};${ga-network_ldflags}")
            target_compile_options(ga::ga INTERFACE ${ga-cflags} ${ga-cppflags})
        else (ga-config)
            message("Local GA not found, attempting to fetch")
            include(FetchContent)
            FetchContent_Declare(
                    ga
                    GIT_REPOSITORY https://github.com/GlobalArrays/ga.git
                    GIT_TAG 2518e23433385bfa3726d507b8cd0d7ed038021b
            )
            FetchContent_MakeAvailable(ga)
            foreach (t ga ga_src dra eaf elio sf linalg)
                install(TARGETS ${t} EXPORT gaTargets
                        LIBRARY DESTINATION lib
                        ARCHIVE DESTINATION lib
                        INCLUDES DESTINATION include
                        PUBLIC_HEADER DESTINATION include
                        )
            endforeach ()
            install(EXPORT gaTargets
                    FILE gaTargets.cmake
                    NAMESPACE ga::
                    DESTINATION lib/cmake/ga)
        endif (ga-config)
    endif ()
endif ()
