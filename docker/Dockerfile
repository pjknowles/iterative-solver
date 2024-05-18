FROM python
#FROM ubuntu
#RUN apt update
#RUN DEBIAN_FRONTEND=noninteractive apt install -y cmake git g++ gfortran doxygen graphviz bash rsync curl mpich libblas-dev liblapack-dev liblapacke-dev libeigen3-dev libhdf5-dev libhdf5-mpich-dev clang ccache make python3 python-is-python3 python3-pip python3-build python3-venv wget ninja-build python3-virtualenv
#RUN DEBIAN_FRONTEND=noninteractive apt install -y python3
#RUN apt update
#RUN DEBIAN_FRONTEND=noninteractive apt install -y python3 python-is-python3 python3-pip python3-build python3-venv
#ENV CMAKE_GENERATOR=Ninja
#RUN export VERSION=5.8.2 && wget https://github.com/GlobalArrays/ga/releases/download/v${VERSION}/ga-${VERSION}.tar.gz && tar xzf ga-${VERSION}.tar.gz && cd ga-${VERSION} && CFLAGS=-fpic ./configure --disable-f77 --with-mpi3 --without-blas --without-lapack && make install && cd .. && rm -rf ga-${VERSION}*

RUN pip install build
RUN pip install pylib-openblas
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y libopenmpi-dev libopenblas-dev liblapacke-dev libeigen3-dev
#RUN pip install openmpi