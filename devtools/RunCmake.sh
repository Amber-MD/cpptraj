#!/bin/bash

# Can be used to test the cmake install.
# Assumes 'git submodule update --init --recursive' has been run in top dir.

if [ -z "$CPPTRAJHOME" ] ; then
  echo "CPPTRAJHOME must be set."
  exit 1
fi
HOME=$CPPTRAJHOME

installdir=$HOME/install
if [ ! -d "$installdir" ] ; then
  mkdir $installdir
fi

#export BUILD_FLAGS="-DOPENMP=TRUE"
#export BUILD_FLAGS="-DMPI=TRUE ${BUILD_FLAGS}"
#-DNetCDF_LIBRARIES_C=$HOME/lib/libnetcdf.so -DNetCDF_INCLUDES=$HOME/include
COMPILER=gnu
cmake .. $BUILD_FLAGS -DCOMPILER=${COMPILER^^} -DINSTALL_HEADERS=FALSE \
         -DCMAKE_INSTALL_PREFIX=$installdir -DCMAKE_LIBRARY_PATH=$HOME/lib \
         -DPRINT_PACKAGING_REPORT=TRUE 

