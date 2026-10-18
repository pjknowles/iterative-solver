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
# Every compiler and library comes from conda. apt supplies tools only, so nothing is
# resolved against a second package manager: when apt's compilers, MPI, BLAS, Eigen
# and HDF5 were installed alongside conda's, the build used conda's compilers with
# whichever libraries CMake happened to find, and Global Arrays was compiled with
# apt's gcc against apt's MPICH but then run against conda's MPICH. Its counterpart
# Dockerfile.debian is all-apt.
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    bash ca-certificates ccache curl doxygen git graphviz make ninja-build openssh-client rsync \
    valgrind wget \
 && rm -rf /var/lib/apt/lists/*
ENV CMAKE_GENERATOR=Ninja
COPY python/requirements.txt .
# The base image already ships conda-libmamba-solver as its default solver, so
# nothing is installed to get it. Reinstalling it with the classic solver, as this
# did, ran for over 45 minutes and was then killed for lack of memory.
RUN CONDA_PLUGINS_AUTO_ACCEPT_TOS=true conda install --file requirements.txt -y -c conda-forge mkl mkl-include hdf5 gcc=14 fortran-compiler
# Global Arrays, built last, with conda's mpicc against conda's MPICH: the toolchain and
# the MPI that the project is then built and tested with. Its configure insists on finding
# a remote shell, which is why openssh-client is installed above.
RUN export VERSION=5.9.2 \
 && wget -q https://github.com/GlobalArrays/ga/releases/download/v${VERSION}/ga-${VERSION}.tar.gz \
 && tar xzf ga-${VERSION}.tar.gz \
 && cd ga-${VERSION} \
 && ./configure --disable-f77 --with-mpi3 --without-blas --without-lapack MPICC=mpicc \
      CFLAGS="-O2 -Wno-error=incompatible-pointer-types" \
 && make install \
 && cd .. && rm -rf ga-${VERSION}*
