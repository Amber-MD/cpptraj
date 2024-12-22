#!/bin/sh

arch=`uname`

if [ "$arch" == "Darwin" ]; then
   ./configure -openblas -arpack -fftw3 -zlib -bzlib -shared clang
   make install
else
   ./configure -readline -openmm -openblas -fftw3 -zlib -bzlib -shared -arpack gnu
   make install
   make clean
   ./configure -readline -openmm -openblas -fftw3 -zlib -bzlib -shared -mpi -arpack gnu
   make install
fi


rsync -a README.md config.h LICENSE dat bin lib test $PREFIX

# mkdir -p $PREFIX/doc
# rsync -a doc/cpptraj.pdf $PREFIX/doc

# Export CPPTRAJHOME automatically
mkdir -p ${PREFIX}/etc/conda/{activate,deactivate}.d
cp ${RECIPE_DIR}/activate.sh ${PREFIX}/etc/conda/activate.d/cpptraj.sh
cp ${RECIPE_DIR}/deactivate.sh ${PREFIX}/etc/conda/deactivate.d/cpptraj.sh
