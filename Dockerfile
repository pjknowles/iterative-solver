FROM ubuntu:impish-20210827
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake git g++ gfortran doxygen graphviz bash rsync curl mpich libblas-dev liblapack-dev liblapacke-dev libeigen3-dev libhdf5-dev libhdf5-mpich-dev clang
RUN git clone https://github.com/GlobalArrays/ga && cd ga && git checkout 2518e23433385bfa3726d507b8cd0d7ed038021b && ./autogen.sh && ./configure --disable-f77 --with-mpi3 && make install
