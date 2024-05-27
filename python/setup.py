import os

import numpy
from setuptools import setup, Extension
from Cython.Build import cythonize

import platform

if platform.system() == "Darwin":
    extra_args = ['-std=c++17', "-mmacosx-version-min=13"]
elif platform.system() == "Windows":
    extra_args = ['/std:c++17']
else:
    extra_args = ['-std=c++17']
prefix = os.environ['PREFIX'] if 'PREFIX' in os.environ else os.environ['CONDA_PREFIX']
iterative_solver_prefix = os.environ['ITERATIVE_SOLVER_PREFIX'] if 'ITERATIVE_SOLVER_PREFIX' in os.environ else \
os.environ['PREFIX']

ext = Extension('iterative_solver.iterative_solver_extension',
                sources=["iterative_solver/iterative_solver_extension.pyx"],
                language="c++",
                extra_compile_args=extra_args,
                include_dirs=[numpy.get_include(), iterative_solver_prefix + '/include', prefix + '/include'],
                define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
                libraries=["blas", "mpi"],
                extra_objects=[iterative_solver_prefix + "/lib/lib" + n + ".a" for n in
                               ["iterative-solver", "profiler", "utilities"]]
                )

setup(
    name="iterative_solver",
    version=(os.environ['VERSION']),
    license="MIT",
    packages=['iterative_solver'],
    ext_modules=cythonize(
        [ext],
        language_level=3,
    ),
)
