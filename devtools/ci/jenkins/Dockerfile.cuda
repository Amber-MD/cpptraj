FROM nvidia/cuda:10.0-devel-ubuntu18.04

ENV CUDA_HOME=/usr/local/cuda

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y gfortran gcc g++ bzip2 libfftw3-dev automake make libbz2-dev openmpi-bin \
                       openmpi-common libopenmpi-dev zlib1g-dev netcdf-bin libnetcdf-dev \
                       liblapack-dev libarpack2-dev