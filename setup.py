# import os

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import subprocess
import shutil
import pathlib

root_dir_ = pathlib.Path(__file__).parent.resolve()
python_source_dir_ = root_dir_ / 'src' / 'iterative_solver'
build_dir_ = root_dir_ / 'build'
cmake_build_dir_ = build_dir_ / 'cmake-build'
cmake_build_dir_.mkdir(parents=True, exist_ok=True)

subprocess.run(
    [
        'cmake',
        '-DCMAKE_CXX_FLAGS=-fPIC',
        '-DCMAKE_BUILD_TYPE=Release',
        '-DDEPENDENCYMANAGER_FETCHCONTENT=OFF', '-DLINEARALGEBRA_ARRAY_HDF5=OFF', '-DLINEARALGEBRA_ARRAY_GA=OFF',
        '-DFORTRAN=OFF',
        '-DBUILD_SHARED_LIBS=OFF', '-DLIBRARY_ONLY=ON',
        # '-DCMAKE_MACOSX_RPATH=OFF',
        # '--trace',
        '-S', str(root_dir_), '-B', str(cmake_build_dir_),
    ],
    # env=dict(os.environ)|{
    #     'CXX':'/Users/peterk/.local/mpich/bin/mpicxx',
    #     'CC':'/Users/peterk/.local/mpich/bin/mpicc',
    #      },
    shell=False)
subprocess.run(['cmake', '--build', str(cmake_build_dir_), '-v', '--config', 'Release'], shell=False)
for suffix in ['so', 'dylib', 'a']:
    for path in pathlib.Path(cmake_build_dir_).glob('**/*.' + suffix):
        shutil.copy2(str(path), str(cmake_build_dir_.parent))

ext = Extension('iterative_solver',
                sources=["src/iterative_solver/iterative_solver.pyx"],
                language="c++",
                include_dirs=[numpy.get_include(), str(root_dir_ / 'src'), str(root_dir_ / 'dependencies/profiler/src'),
                              str(root_dir_ / 'dependencies/utilities/src'),
                              str(python_source_dir_),
                              ],
                extra_compile_args=["-std=c++17"],
                extra_objects=[
                    "build/libiterative-solver.a",
                    "build/libutilities.a",
                    "build/libprofiler.a",
                    # "/Users/peterk/.local/mpich/lib/libmpi.a",
                ],
                define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
                libraries=["iterative-solver", "utilities", "blas", "mpi"],
                library_dirs=[str(cmake_build_dir_.parent)]
                )

setup(
    name="iterative_solver",
    ext_modules=cythonize(
        [ext],
        language_level=3,
    ),
)
