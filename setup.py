import os

from setuptools import setup, Extension
# from distutils.core import Extension
# from Cython.Distutils.extension import Extension
# from distutils.core import setup
from Cython.Build import cythonize
import numpy
import subprocess
import shutil
import pathlib

root_dir_ = pathlib.Path(__file__).parent.resolve()
python_source_dir_ = root_dir_ / 'python' / 'molpro' / 'iterative_solver'
python_source_dir_ = root_dir_ / 'src' / 'iterative_solver'
build_dir_ = root_dir_ / 'build'
cmake_build_dir_ = build_dir_ / 'cmake-build'
cmake_build_dir_.mkdir(parents=True, exist_ok=True)
try:
    # os.remove(cmake_build_dir_ / 'CMakeCache.txt')
    # shutil.rmtree(cmake_build_dir_)
    pass
except:
    pass
# print('root_dir_:', root_dir_)
# subprocess.run(['ls','-lR',root_dir_/'dependencies'])
# subprocess.run(['tar','cf','/tmp/iterative-solver.tar','-C',root_dir_.parent,'.'])
# print('cmake_build_dir_:', cmake_build_dir_)
subprocess.run(
    [
        'cmake',
        '-DCMAKE_BUILD_TYPE=Release',
        '-DDEPENDENCYMANAGER_FETCHCONTENT=OFF', '-DLINEARALGEBRA_ARRAY_HDF5=OFF', '-DFORTRAN=OFF',
        '-DBUILD_SHARED_LIBS=ON',
        # '-DCMAKE_MACOSX_RPATH=OFF',
        # '--trace',
        '-S', str(root_dir_), '-B', str(cmake_build_dir_),
    ],
    shell=False)
# subprocess.run(['ls','-lR',root_dir_/'dependencies'])
# subprocess.run(['find',cmake_build_dir_,'-name','LibraryManager.cmake','-ls'])
# subprocess.run(['grep','-ir','library-manager',cmake_build_dir_])
subprocess.run(['cmake', '--build', str(cmake_build_dir_), '-v', '--config', 'Release'], shell=False)
for suffix in ['so', 'dylib']:
    for path in pathlib.Path(cmake_build_dir_).glob('**/*.' + suffix):
        shutil.copy2(str(path), str(cmake_build_dir_.parent))

print('after cmake build')
ext = Extension('iterative_solver',
                # sources=["python/molpro/iterative_solver/iterative_solver.pyx","python/molpro/iterative_solver/iterative_solver.pxd"],
                sources=["src/iterative_solver/iterative_solver.pyx"],
                # sources=["iterative_solver.pyx"],
                language="c++",
                include_dirs=[numpy.get_include(), str(root_dir_ / 'src'), str(root_dir_ / 'dependencies/profiler/src'),
                              str(root_dir_ / 'dependencies/utilities/src'),
                              str(python_source_dir_),
                              # str(python_source_dir_ / 'src' / 'molpro'/'iterative_solver'),
                              ],
                extra_compile_args=["-std=c++17"],
                define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
                libraries=["iterative-solver", "utilities"],
                library_dirs=[str(cmake_build_dir_.parent)]
                )

print('before setup',os.getcwd())
subprocess.run(['ls','-lR'])
setup(
    name="iterative_solver",
    ext_modules=cythonize(
        [ext],
        language_level=3,
    ),
    #    package_dir = {
    #        "iterative_solver" : "python/molpro/iterative_solver"
    #    }
)
print('after setup')
