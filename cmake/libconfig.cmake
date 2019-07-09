include(FetchContent)

# Declare an external git-hosted library on which this project depends
# The first parameter is a character string that will be used to generate file
# and target names, and as a handle for a subsequent get_dependency() call.
# The second parameter is the URL where the git repository can be found.
# The node in the git repository is specified separately in the file
# ${CMAKE_SOURCE_DIR}/dependencies/${NAME}_SHA1
function(declare_dependency NAME URL)
    set(_private_dependency_${NAME}_directory "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "dependency directory for ${NAME}")
    get_dependency_name(${NAME})
    file(STRINGS "${_SHA_file}" GIT_TAG)
    message(STATUS "Declare dependency NAME=${NAME} URL=${URL} TAG=${GIT_TAG} DEPENDENCY=${_dependency_name}")
    FetchContent_Declare(
            ${_dependency_name}
            SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/${NAME}"
            GIT_REPOSITORY ${URL}
            GIT_TAG ${GIT_TAG}
            UPDATE_COMMAND "${CMAKE_SOURCE_DIR}/dependencies/checkout.sh" ${NAME} ${CMAKE_BUILD_TYPE}
    )
endfunction()

# Load an external git-hosted library on which this project depends.
# The first parameter is the name of the dependency, as passed previously to
# declare_dependency().
# A second parameter can be given, which if evaluating to true includes the
# library in the cmake build via add_subdirectory(). If omitted, true is assumed.
function(get_dependency name)
    get_dependency_name(${name})
    if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.15)
        message(VERBOSE "get_dependency(${name})")
    endif ()
    FetchContent_GetProperties(${_dependency_name})
    if (NOT ${_dependency_name}_POPULATED)
        file(LOCK "${CMAKE_SOURCE_DIR}/dependencies/.${name}_lockfile" GUARD FILE TIMEOUT 1000)
        FetchContent_Populate(${_dependency_name})
        file(LOCK "${CMAKE_SOURCE_DIR}/dependencies/.${name}_lockfile" RELEASE)
    endif ()
    if (ARGV1 OR (NOT DEFINED ARGV1))
        if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.15)
            message(VERBOSE "add_subdirectory() for get_dependency ${name}")
        endif ()
        add_subdirectory(${${_dependency_name}_SOURCE_DIR} ${${_dependency_name}_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif ()
    foreach (s SOURCE_DIR BINARY_DIR POPULATED)
        set(${_dependency_name}_${s} "${${_dependency_name}_${s}}" PARENT_SCOPE)
    endforeach ()
endfunction()

# Completion of configuration for a library, including installation and package
# configuration. After installation, ${CMAKE_INSTALL_PREFIX}/bin/${LIBRARY_NAME}-config
# is a script that can be used to discover build options for using the library.
# LIBRARY_NAME: target name of library, which should already have been defined
# via add_library()
# DEPENDENCIES: A (possibly empty) list of existing target names of libraries on
# which LIBRARY_NAME depends. For each ${lib} in the list, the variable
# ${DEPENDENCY_${lib}} contains additional specification to be appended to the
# library name in the export definition expressed through the find_dependency()
# command for constructing the package configuration file. Normally, this will
# be the version number of the dependency.
function(configure_library LIBRARY_NAME DEPENDENCIES)
    if (DEPENDENCIES)
        message(STATUS "Configure_library ${LIBRARY_NAME} with dependencies ${DEPENDENCIES}")
    else ()
        message(STATUS "Configure_library ${LIBRARY_NAME}")
    endif ()
    string(TOUPPER ${LIBRARY_NAME} PROJECT_UPPER_NAME)
    add_library(${LIBRARY_NAME}::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})
    target_include_directories(${LIBRARY_NAME} PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
            $<INSTALL_INTERFACE:include>
            )
    target_compile_definitions(${LIBRARY_NAME} PRIVATE NOMAIN)
    if (Molpro_SOURCE_DIR)
        set(MOLPRO 1)
        target_include_directories(${LIBRARY_NAME} PRIVATE "${CMAKE_BINARY_DIR}/src" "${Molpro_SOURCE_DIR}/build" "${Molpro_SOURCE_DIR}/src")
    endif ()
    if (FORTRAN)
        set(${PROJECT_UPPER_NAME}_FORTRAN 1)
        if (INTEGER8)
            set(${PROJECT_UPPER_NAME}_I8 1)
        endif ()
        include(CheckFortranCompilerFlag)
        CHECK_Fortran_COMPILER_FLAG("-Wall" _Wallf)
        if (_Wallf)
            set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall")
        endif ()
    endif ()

    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-Wall" _Wall)
    if (_Wall)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    endif ()
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}-config.h.in" "
#ifndef ${PROJECT_UPPER_NAME}_CONFIG_H
#define ${PROJECT_UPPER_NAME}_CONFIG_H
#cmakedefine ${PROJECT_UPPER_NAME}_FORTRAN
#cmakedefine ${PROJECT_UPPER_NAME}_I8
#ifndef MOLPRO
#cmakedefine MOLPRO
#endif
#endif
")
    configure_file("${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}-config.h.in" ${LIBRARY_NAME}-config.h)

    install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

    install(TARGETS ${LIBRARY_NAME} EXPORT ${LIBRARY_NAME}Targets LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            RUNTIME DESTINATION bin
            INCLUDES DESTINATION include
            PUBLIC_HEADER DESTINATION include
            )
    foreach (dep ${DEPENDENCIES})
        install(TARGETS ${dep} EXPORT ${LIBRARY_NAME}Targets LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
                RUNTIME DESTINATION bin
                INCLUDES DESTINATION include
                PUBLIC_HEADER DESTINATION include
                )
    endforeach ()
    install(EXPORT ${LIBRARY_NAME}Targets
            FILE ${LIBRARY_NAME}Targets.cmake
            NAMESPACE ${LIBRARY_NAME}::
            DESTINATION lib/cmake/${LIBRARY_NAME}
            )

    include(CMakePackageConfigHelpers)
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}Config.cmake"
            "include(CMakeFindDependencyMacro)
")
    foreach (dep ${DEPENDENCIES})
        file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}Config.cmake"
                "find_dependency(${dep} ${DEPENDENCY_${dep}})
")
    endforeach ()
    file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}Config.cmake" "
include(\"\${CMAKE_CURRENT_LIST_DIR}/${LIBRARY_NAME}Targets.cmake\")
")
    write_basic_package_version_file("${LIBRARY_NAME}ConfigVersion.cmake"
            VERSION ${CMAKE_PROJECT_VERSION}
            COMPATIBILITY SameMajorVersion
            )
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}Config.cmake" "${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}ConfigVersion.cmake"
            DESTINATION lib/cmake/${LIBRARY_NAME}
            )

    set(CONFIG_CPPFLAGS "-I${CMAKE_INSTALL_PREFIX}/include")
    get_target_property(FLAGS ${LIBRARY_NAME} INTERFACE_COMPILE_DEFINITIONS)
    if (FLAGS)
        foreach (flag ${FLAGS})
            set(CONFIG_CPPFLAGS "${CONFIG_CPPFLAGS} -D${flag}")
        endforeach ()
    endif ()
    #set(CONFIG_FCFLAGS "${CMAKE_Fortran_MODDIR_FLAG}${CMAKE_INSTALL_PREFIX}/include")
    set(CONFIG_FCFLAGS "-I${CMAKE_INSTALL_PREFIX}/include") #TODO should not be hard-wired -I
    set(CONFIG_LDFLAGS "-L${CMAKE_INSTALL_PREFIX}/lib")
    set(CONFIG_LIBS "-l${LIBRARY_NAME}")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}/config.in" "
#!/bin/sh

prefix=\"@CMAKE_INSTALL_PREFIX@\"
exec_prefix=\"\${prefix}\"

usage='Usage: '$0' [--cppflags] [--cxxflags] [--exec-prefix] [--ldflags] [--libs] [--prefix] [--version]

 For example, 'test.cpp' may be compiled to produce 'test' as follows:

  c++ -o test test.cpp `'$0' --cppflags --cxxflags --ldflags --libs`'

if test $# -eq 0; then
      echo \"\${usage}\" 1>&2
      exit 1
fi

while test $# -gt 0; do
  case \"\$1\" in
    -*=*) optarg=`echo \"\$1\" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
    *) optarg= ;;
  esac
  case $1 in
    --prefix=*)
      prefix=$optarg
      ;;
    --prefix)
      echo $prefix
      ;;
    --exec-prefix=*)
      exec_prefix=$optarg
      ;;
    --exec-prefix)
      echo $exec_prefix
      ;;
    --version)
      echo '@CMAKE_PROJECT_VERSION@'
      ;;
    --cflags)
      echo '@CONFIG_CFLAGS@'
      ;;
    --cxxflags)
      echo '@CONFIG_CXXFLAGS@'
      ;;
    --cppflags)
      echo '@CONFIG_CPPFLAGS@'
      ;;
    --fflags)
      echo '@CONFIG_FFLAGS@'
      ;;
    --fcflags)
      echo '@CONFIG_FCFLAGS@'
      ;;
    --ldflags)
      echo '@CONFIG_LDFLAGS@'
      ;;
    --libs)
      echo '@CONFIG_LIBS@'
      ;;
    *)
      echo \"\${usage}\" 1>&2
      exit 1
      ;;
  esac
  shift
done


")
    configure_file("${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}/config.in" ${LIBRARY_NAME}-config @ONLY)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}-config DESTINATION bin)
endfunction()

function(get_dependency_name dep)
    set(_dependency_name _private_dep_${dep} PARENT_SCOPE)
    set(_SHA_file "${_private_dependency_${NAME}_directory}/${dep}_SHA1" PARENT_SCOPE)
endfunction()

