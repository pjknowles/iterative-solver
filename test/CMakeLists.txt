include(FetchContent)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.10.0"
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)



add_subdirectory(array)
add_subdirectory(itsolv)
