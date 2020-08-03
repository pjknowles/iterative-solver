# parent image
FROM ubuntu:20.04

# Install any needed packages
RUN apt-get update 
RUN apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake git g++ gfortran doxygen graphviz bash rsync curl
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake mpich
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake libblas-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake liblapack-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake libeigen3-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake libhdf5-dev libhdf5-mpich-dev
RUN git clone https://github.com/GlobalArrays/ga
RUN cd ga; git checkout 2518e23433385bfa3726d507b8cd0d7ed038021b; ./autogen.sh
RUN cd ga && ./configure --disable-f77 --with-mpi3 && make install
