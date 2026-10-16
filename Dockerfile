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
RUN CONDA_PLUGINS_AUTO_ACCEPT_TOS=true conda install --solver=classic conda-forge::conda-libmamba-solver conda-forge::libmamba conda-forge::libmambapy conda-forge::libarchive
RUN CONDA_PLUGINS_AUTO_ACCEPT_TOS=true conda install --file requirements.txt -y -c conda-forge mkl mkl-include gcc=14 fortran-compiler
