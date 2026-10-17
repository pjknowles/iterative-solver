# Pinned by digest so that a CI result only changes when this repository does.
# This is the image that docker.io/continuumio/miniconda3:latest resolved to on
# 2026-10-16 (built 2026-08-28), so pinning it changes nothing about what CI runs
# today. Bump it deliberately; note that the newest version tag (v25.11.1-1) is a
# different image from latest.
FROM continuumio/miniconda3@sha256:eca594d684f495c1a02beff33a9fab53aec8c5830eaf431bb149912dc6c9e4c1
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
