# parent image
FROM ubuntu:groovy-20200812

# Install any needed packages
RUN apt-get update 
RUN apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake git g++ gfortran doxygen graphviz bash rsync curl
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y mpich
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libblas-dev liblapack-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libeigen3-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libhdf5-dev libhdf5-mpich-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y clang
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y liblapacke-dev
RUN git clone https://github.com/GlobalArrays/ga && cd ga && git checkout 2518e23433385bfa3726d507b8cd0d7ed038021b && ./autogen.sh && ./configure --disable-f77 --with-mpi3 && make install
