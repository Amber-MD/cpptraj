FROM ubuntu:18.04

RUN mkdir -p /opt/intel && \
    mkdir -p /opt/pgi && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y gfortran gcc g++ bzip2 libfftw3-dev automake make libbz2-dev \
                       mpich libmpich-dev zlib1g-dev netcdf-bin libnetcdf-dev \
                       liblapack-dev libarpack2-dev

# Set up for the Intel compilers. This won't do anything if they're not mounted
ENV PATH=/opt/intel/compilers_and_libraries/linux/bin/intel64:/opt/intel/openmpi/bin:${PATH}
ENV MKL_HOME=/opt/intel/compilers_and_libraries/linux/mkl
ENV LD_LIBRARY_PATH=${MKL_HOME}/lib/intel64:/opt/intel/compilers_and_libraries/linux/lib/intel64:/opt/intel/openmpi/lib

# Set up for the PGI compilers. This won't do anything if they're not mounted
ENV PATH=/opt/pgi/linux86-64/2018/bin:/opt/pgi/linux86-64/2018/mpi/bin:${PATH}
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/pgi/linux86-64/18.10/lib:/opt/pgi/linux86-64/2018/mpi/lib
