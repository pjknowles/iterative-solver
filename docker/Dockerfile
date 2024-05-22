#FROM ubuntu
FROM continuumio/miniconda3
RUN apt update
#RUN DEBIAN_FRONTEND=noninteractive apt install -y cmake git g++ gfortran doxygen graphviz bash rsync curl mpich libblas-dev liblapack-dev liblapacke-dev libeigen3-dev libhdf5-dev libhdf5-mpich-dev clang ccache make python3 python-is-python3 python3-pip python3-build python3-venv wget ninja-build python3-virtualenv
RUN DEBIAN_FRONTEND=noninteractive apt install -y libeigen3-dev python3 python-is-python3 python3-pip python3-build python3-venv wget ninja-build python3-virtualenv

#RUN  wget -nv https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-$(arch).sh -O /tmp/miniconda-installer.sh 
#RUN bash /tmp/miniconda-installer.sh -b -p /opt/conda
#RUN  . /opt/conda/etc/profile.d/conda.sh

SHELL ["/bin/bash", "--login", "-c"]
RUN conda init bash
RUN conda create -c conda-forge -n iterative-solver -y c-compiler cxx-compiler openblas cython cmake
RUN conda install -c conda-forge -y -n iterative-solver mpich
RUN conda install -c conda-forge -y -n iterative-solver numpy
RUN conda install -c conda-forge -y -n iterative-solver lapack liblapacke
RUN echo 'conda activate iterative-solver' >> ~/.bashrc

ENV CMAKE_GENERATOR=Ninja
#ENV CONDA_DEFAULT_ENV=iterative-solver
#RUN export VERSION=5.8.2 && wget https://github.com/GlobalArrays/ga/releases/download/v${VERSION}/ga-${VERSION}.tar.gz && tar xzf ga-${VERSION}.tar.gz && cd ga-${VERSION} && CFLAGS=-fpic ./configure --disable-f77 --with-mpi3 --without-blas --without-lapack && make install && cd .. && rm -rf ga-${VERSION}*