FROM continuumio/miniconda:latest

RUN mkdir /cpptraj && \
    mkdir /.conda /.local && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y gfortran gcc g++ bzip2 libfftw3-dev automake make libbz2-dev \
                       mpich libmpich-dev zlib1g-dev netcdf-bin libnetcdf-dev \
                       liblapack-dev libarpack2-dev

# Make it so Jenkins UIDs can install python packages to miniconda
RUN chmod -R a+w /opt/conda && \
    chmod -R a+w /.conda && \
    chmod -R a+w /.local

COPY . /cpptraj

ENV CPPTRAJHOME "/cpptraj"

RUN cd /cpptraj && \
    ./configure -openmp -shared gnu && \
    make -j4 libcpptraj

