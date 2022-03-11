#!/bin/bash

# Can be used to test the cmake install.

WORKDIR=`pwd`
# For safety
if [ "`basename $WORKDIR`" != 'build' ] ; then
  echo "Please run from a 'build' directory"
  exit 1
fi

if [ -z "$CPPTRAJHOME" ] ; then
  echo "CPPTRAJHOME must be set."
  exit 1
fi
HOME=$CPPTRAJHOME

installdir=$HOME/install
if [ ! -d "$installdir" ] ; then
  mkdir $installdir
fi

CUDA=0
OPENMP=0
MPI=0
FFT=1
ARPACK=1
COLOR=1
COMPILER=gnu
while [ ! -z "$1" ] ; do
  case "$1" in
    'cuda' ) CUDA=1 ;;
    'openmp' ) OPENMP=1 ;;
    'mpi' ) MPI=1 ;;
    'all' ) CUDA=1 ; OPENMP=1 ; MPI=1 ;;
    'nofftw3') FFT=0 ;;
    'noarpack') ARPACK=0 ;;
    'nocolor') COLOR=0 ;;
    'intel' ) COMPILER=intel ;;
    'autotest' )
      FFT=0
      ARPACK=0
      COMPILER='AUTO'
      #export BUILD_FLAGS="-DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc ${BUILD_FLAGS}"
      #export BUILD_FLAGS="-DCMAKE_Fortran_COMPILER=ifort ${BUILD_FLAGS}"
    ;;
  esac
  shift
done

if [ $CUDA -eq 1 ] ; then
  export BUILD_FLAGS="-DCUDA=TRUE ${BUILD_FLAGS}"
fi
if [ $OPENMP -eq 1 ] ; then
  export BUILD_FLAGS="-DOPENMP=TRUE ${BUILD_FLAGS}"
fi
if [ $MPI -eq 1 ] ; then
  export BUILD_FLAGS="-DMPI=TRUE ${BUILD_FLAGS}"
  export BUILD_FLAGS="-DPnetCDF_LIBRARY=$PNETCDF_HOME/lib/libpnetcdf.a ${BUILD_FLAGS}"
  export BUILD_FLAGS="-DPnetCDF_INCLUDE_DIR=$PNETCDF_HOME/include ${BUILD_FLAGS}"
fi
if [ $FFT -eq 0 ] ; then
  export BUILD_FLAGS="-DUSE_FFT=FALSE ${BUILD_FLAGS}"
fi
if [ $ARPACK -eq 0 ] ; then
  export BUILD_FLAGS="-DUSE_ARPACK=FALSE ${BUILD_FLAGS}"
fi
if [ $COLOR -eq 0 ] ; then
  export BUILD_FLAGS="-DCOLOR_CMAKE_MESSAGES=FALSE ${BUILD_FLAGS}"
fi
#-DNetCDF_LIBRARIES_C=$HOME/lib/libnetcdf.so -DNetCDF_INCLUDES=$HOME/include
cmake .. \
         $BUILD_FLAGS -DCOMPILER=${COMPILER^^} -DINSTALL_HEADERS=FALSE \
         -DCMAKE_INSTALL_PREFIX=$installdir -DCMAKE_LIBRARY_PATH=$HOME/lib \
         -DPRINT_PACKAGING_REPORT=TRUE 

