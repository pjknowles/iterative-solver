from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import subprocess
import shutil
import pathlib

current_ = pathlib.Path(__file__).parent.resolve()
build_ = (current_ / 'build' / 'cmake-build')
build_.mkdir(parents=True, exist_ok=True)
source_ = current_.parent.parent.parent.resolve()
subprocess.run(
    [
        'cmake',
        '-DCMAKE_BUILD_TYPE=Release',
        '-DDEPENDENCYMANAGER_FETCHCONTENT=ON', '-DLINEARALGEBRA_ARRAY_HDF5=OFF', '-DFORTRAN=OFF',
        '-DBUILD_SHARED_LIBS=ON',
        '-S', str(source_), '-B', str(build_),
    ],
    shell=False)
subprocess.run(['cmake', '--build', str(build_), '--config', 'Release'], shell=False)
for suffix in ['so', 'dylib']:
    for path in pathlib.Path(build_).glob('**/*.'+suffix):
        shutil.copy2(str(path), str(build_.parent))

ext = Extension('iterative_solver',
                sources=["iterative_solver.pyx"],
                language="c++",
                include_dirs=[numpy.get_include(), '../../../src', '../../../dependencies/profiler/src',
                              '../../../dependencies/utilities/src'],
                extra_compile_args=["-std=c++17"],
                define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
                libraries=["iterative-solver", "utilities"],
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
