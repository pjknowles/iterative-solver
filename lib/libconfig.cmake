string(TOUPPER ${PROJECT_NAME} PROJECT_UPPER_NAME)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
        $<INSTALL_INTERFACE:include>
        )
target_compile_definitions(${PROJECT_NAME} PRIVATE NOMAIN)
if (Molpro_SOURCE_DIR)
    set(MOLPRO 1)
    target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_BINARY_DIR}" "${CMAKE_BINARY_DIR}/src" "${Molpro_SOURCE_DIR}/build" "${Molpro_SOURCE_DIR}/src")
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
configure_file(${PROJECT_NAME}-config.h.in ${PROJECT_NAME}-config.h)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        PUBLIC_HEADER DESTINATION include
        )
install(EXPORT ${PROJECT_NAME}Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION lib/cmake/${PROJECT_NAME}
        )

include(CMakePackageConfigHelpers)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        "include(CMakeFindDependencyMacro)
")
foreach (dep ${DEPENDENCIES})
    file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
            "find_dependency(${dep} ${DEPENDENCY_${dep}})
")
endforeach ()
file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" "
include(\"\${CMAKE_CURRENT_LIST_DIR}/${PROJECT_NAME}Targets.cmake\")
")
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
        VERSION ${CMAKE_PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
        )
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION lib/cmake/${PROJECT_NAME}
        )

set(CONFIG_CPPFLAGS "-I${CMAKE_INSTALL_PREFIX}/include")
get_target_property(FLAGS ${PROJECT_NAME} INTERFACE_COMPILE_DEFINITIONS)
if (FLAGS)
    foreach (flag ${FLAGS})
        set(CONFIG_CPPFLAGS "${CONFIG_CPPFLAGS} -D${flag}")
    endforeach ()
endif ()
#set(CONFIG_FCFLAGS "${CMAKE_Fortran_MODDIR_FLAG}${CMAKE_INSTALL_PREFIX}/include")
set(CONFIG_FCFLAGS "-I${CMAKE_INSTALL_PREFIX}/include") #TODO should not be hard-wired -I
set(CONFIG_LDFLAGS "-L${CMAKE_INSTALL_PREFIX}/lib")
set(CONFIG_LIBS "-l${PROJECT_NAME}")
configure_file(config.in ${PROJECT_NAME}-config @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config DESTINATION bin)
