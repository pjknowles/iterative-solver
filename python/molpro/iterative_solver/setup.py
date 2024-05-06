# from distutils.core import setup
# from distutils.extension import Extension
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

import subprocess
subprocess.run("""
set -e
mkdir -p build/cmake-build
cd build/cmake-build
cmake -DCMAKE_BUILD_TYPE=Release -DDEPENDENCYMANAGER_FETCHCONTENT=ON -DLINEARALGEBRA_ARRAY_HDF5=OFF -DFORTRAN=OFF -DBUILD_SHARED_LIBS=ON ../../../../..
cmake --build . --config Release
pwd
find . \( -name '*.dylib' -o -name '*.so' \) -exec ls -la {} \;
find . \( -name '*.dylib' -o -name '*.so' \) -exec cp -p {} .. \;
""", shell=True, check=True)

ext = Extension('iterative_solver',
                sources=["iterative_solver.pyx"],
                language="c++",
                include_dirs=[numpy.get_include(), '../../../src','../../../dependencies/profiler/src','../../../dependencies/utilities/src'],
                extra_compile_args=["-std=c++17"],
                define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
                libraries=["iterative-solver","utilities"],
                library_dirs=[
                    "build",
                #     "../../../cmake-build-release-shared/src",
                #     "../../../cmake-build-release-shared/_deps/utilities-build/src",
                ]
                )

setup(
    name="iterative_solver",
    ext_modules=cythonize(
        [ext],
        language_level=3,
    ),
)
