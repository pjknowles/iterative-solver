#ifndef MOLPRO_DEPENDENCIES_LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_FS_H_
#define MOLPRO_DEPENDENCIES_LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_FS_H_
#if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || (defined(__cplusplus) && __cplusplus >= 201703L)) &&            \
    defined(__has_include)

#if __has_include(<filesystem>) && (!defined(__MAC_OS_X_VERSION_MIN_REQUIRED) || __MAC_OS_X_VERSION_MIN_REQUIRED >= 101500)          \
&& ((!defined(__GNUC__) || defined(__llvm__) || defined(__INTEL_COMPILER) || __GNUC__ >= 9))
#define GHC_USE_STD_FS
#include <filesystem>
namespace fs = std::filesystem;
#endif
#endif
#ifndef GHC_USE_STD_FS
#include "../ghc/filesystem.h"
namespace fs = ghc::filesystem;
#endif

#endif // MOLPRO_DEPENDENCIES_LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_FS_H_
