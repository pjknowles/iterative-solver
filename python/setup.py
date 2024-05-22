import os
import platform
import re

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import subprocess
import pathlib

python_dir_ = pathlib.Path(__file__).parent.resolve()
root_dir_ = python_dir_.parent.resolve()
python_source_dir_ = python_dir_ / 'iterative_solver_extension'
build_dir_ = python_dir_ / 'build'
cmake_build_dir_ = python_dir_ / ('cmake-build-' + platform.system() + '-' + platform.machine())
cmake_build_dir_.mkdir(parents=True, exist_ok=True)

subprocess.run(
    [
        'cmake',
        '-DCMAKE_INSTALL_PREFIX=' + os.environ['CONDA_PREFIX'],
        '-DCMAKE_CXX_FLAGS=-fPIC',  # TODO make portable
        '-DCMAKE_BUILD_TYPE=Release',
        '-DDEPENDENCYMANAGER_FETCHCONTENT=OFF',
        '-DLINEARALGEBRA_ARRAY_HDF5=OFF', '-DLINEARALGEBRA_ARRAY_GA=OFF',
        '-DFORTRAN=OFF',
        '-DBUILD_SHARED_LIBS=ON', '-DLIBRARY_ONLY=ON',
        '-S', str(root_dir_), '-B', str(cmake_build_dir_),
    ],
    shell=False)
version_ = '0.0.0'
with open(cmake_build_dir_ / 'project_version.sh', 'r') as f:
    for line in f.readlines():
        if line.strip().startswith('PROJECT_VERSION_FULL'):
            version_ = re.sub('.*=', '', line).strip()
with open('iterative_solver/_version.py', 'w') as f:
    f.write('__version__ = "{}"\n'.format(version_))
subprocess.run(['cmake', '--build', str(cmake_build_dir_), '-t', 'install', '-v', '--config', 'Release'], shell=False)

ext = Extension('iterative_solver_extension',
                sources=["iterative_solver_extension/iterative_solver_extension.pyx"],
                language="c++",
                include_dirs=[numpy.get_include(),
                              str(pathlib.Path(os.environ['CONDA_PREFIX']) / 'include'),
                              str(python_source_dir_),
                              ],
                extra_compile_args=["-std=c++17"],
                define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
                libraries=["blas", "mpi", "iterative-solver", "utilities", "profiler"],
                )

setup(
    name="iterative_solver_extension",
    ext_modules=cythonize(
        [ext],
        language_level=3,
    ),
)
