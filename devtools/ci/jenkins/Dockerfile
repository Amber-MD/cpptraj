ARG BASEIMAGE=ubuntu:18.04
FROM ${BASEIMAGE}}

RUN mkdir -p /opt/intel && \
    mkdir -p /opt/pgi && \
    apt-get update && \
    apt-get install gfortran gcc g++ bzip2 libfftw3-dev automake make \
            openmpi-bin openmpi-common libopenmpi-dev

# Set up for the Intel compilers. This won't do anything if they're not mounted
ENV PATH=${PATH}:/opt/intel/compilers_and_libraries/linux/bin/intel64
ENV LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries/linux/lib/intel64

# Set up for the PGI compilers. This won't do anything if they're not mounted
ENV PATH=${PATH}:/opt/pgi/linux86-64/18.10/bin
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/pgi/linux86-64/18.10/lib
