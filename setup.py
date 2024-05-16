from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import subprocess
import pathlib

root_dir_ = pathlib.Path(__file__).parent.resolve()
python_source_dir_ = root_dir_ / 'src' / 'iterative_solver'
build_dir_ = root_dir_ / 'build'
cmake_build_dir_ = build_dir_ / 'cmake-build'
cmake_build_dir_.mkdir(parents=True, exist_ok=True)

subprocess.run(
    [
        'cmake',
        '-DCMAKE_CXX_FLAGS=-fPIC',  # TODO make portable
        '-DCMAKE_BUILD_TYPE=Release',
        '-DDEPENDENCYMANAGER_FETCHCONTENT=OFF',
        '-DLINEARALGEBRA_ARRAY_HDF5=OFF', '-DLINEARALGEBRA_ARRAY_GA=OFF',
        '-DFORTRAN=OFF',
        '-DBUILD_SHARED_LIBS=OFF', '-DLIBRARY_ONLY=ON',
        '-S', str(root_dir_), '-B', str(cmake_build_dir_),
    ],
    shell=False)
subprocess.run(['cmake', '--build', str(cmake_build_dir_), '-v', '--config', 'Release'], shell=False)

ext = Extension('iterative_solver',
                sources=["src/iterative_solver/iterative_solver.pyx"],
                language="c++",
                include_dirs=[numpy.get_include(), str(root_dir_ / 'src'), str(root_dir_ / 'dependencies/profiler/src'),
                              str(root_dir_ / 'dependencies/utilities/src'),
                              str(python_source_dir_),
                              ],
                extra_compile_args=["-std=c++17"],
                extra_objects=[
                    "build/cmake-build/src/libiterative-solver.a",
                    "build/cmake-build/_deps/utilities-build/src/libutilities.a",
                    "build/cmake-build/_deps/profiler-build/src/libprofiler.a",
                ],
                define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
                libraries=["blas", "mpi"],
                )

setup(
    name="iterative_solver",
    ext_modules=cythonize(
        [ext],
        language_level=3,
    ),
)
