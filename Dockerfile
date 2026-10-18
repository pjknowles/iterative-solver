# Pinned by digest so that a CI result only changes when this repository does.
# This is docker.io/anaconda/miniconda:26.7.1, the successor to the deprecated
# continuumio/miniconda3. It carries the same Miniconda as the continuumio image
# previously pinned here (conda 26.7.1, Python 3.14.7, Debian 13), so the move changes
# where the image comes from, not what is in it. Bump it deliberately.
FROM anaconda/miniconda@sha256:38939b55a6cdca3b5cb24eb8dd9f7f7bcf49cd4ec9e6d72760f8161b906f9085
# This image installs conda in /opt/miniconda3, not /opt/conda, and unlike the
# continuumio one it does not put conda.sh in /etc/profile.d. Without that a login
# shell has no conda function, and test-python.yml's `conda activate base` fails.
RUN ln -s /opt/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y cmake git g++ gfortran doxygen graphviz bash rsync curl mpich libblas-dev liblapack-dev liblapacke-dev libeigen3-dev libhdf5-dev libhdf5-mpich-dev clang ccache ninja-build wget valgrind
ENV CMAKE_GENERATOR=Ninja
RUN export VERSION=5.9.2 && wget https://github.com/GlobalArrays/ga/releases/download/v${VERSION}/ga-${VERSION}.tar.gz && tar xzf ga-${VERSION}.tar.gz && cd ga-${VERSION} && ./configure --disable-f77 --with-mpi3 --without-blas --without-lapack CFLAGS="-O2 -Wno-error=incompatible-pointer-types" && make install && cd .. && rm -rf ga-${VERSION}*
COPY python/requirements.txt .
# The base image already ships conda-libmamba-solver as its default solver, so
# nothing is installed to get it. Reinstalling it with the classic solver, as this
# did, ran for over 45 minutes and was then killed for lack of memory.
# hdf5 comes from conda rather than only apt: the build uses the conda toolchain and
# conda's cmake, which did not find the apt HDF5, so -DLINEARALGEBRA_ARRAY_HDF5=ON
# could not be covered.
RUN CONDA_PLUGINS_AUTO_ACCEPT_TOS=true conda install --file requirements.txt -y -c conda-forge mkl mkl-include hdf5 gcc=14 fortran-compiler
