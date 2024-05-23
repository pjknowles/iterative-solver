#!/bin/sh
set -euo pipefail

if [ -z "$CONDA_PREFIX" -o $(basename $CONDA_PREFIX) = "base" ]; then
  echo $0: must run this script in a non-base Conda environment
  echo CONDA_PREFIX=$CONDA_PREFIX $(basename $CONDA_PREFIX)
  exit 1
fi

python_dir=$(realpath $(dirname $0))
root_dir=$(realpath $python_dir/..)
conda install -c conda-forge -y --file $python_dir/requirements.txt
cmake_build_dir=$python_dir/cmake-build-$(uname)-$(uname -m)
mkdir -p $cmake_build_dir

cmake \
  -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -DCMAKE_CXX_FLAGS=-fPIC \
  -DCMAKE_BUILD_TYPE=Release \
  -DDEPENDENCYMANAGER_FETCHCONTENT=OFF \
  -DLINEARALGEBRA_ARRAY_HDF5=OFF -DLINEARALGEBRA_ARRAY_GA=OFF \
  -DFORTRAN=OFF \
  -DBUILD_SHARED_LIBS=ON, -DLIBRARY_ONLY=ON \
  -S $root_dir -B $cmake_build_dir

export VERSION=$(grep PROJECT_VERSION= $cmake_build_dir/project_version.sh | sed -e 's/.*=//')
if [ -z "$VERSION" ]; then
  export VERSION='0.0.0'
fi
echo VERSION=$VERSION
echo '__version__ = "'$VERSION'"' >$python_dir/iterative_solver/_version.py
cmake --build $cmake_build_dir -t install -v --config Release

python -m pip install -e $python_dir
