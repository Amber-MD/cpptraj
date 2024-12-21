#!/bin/sh

arch=`uname`

if [ "$arch" == "Darwin" ]; then
   ./configure -openblas -arpack -fftw3 -zlib -bzlib -shared clang
else
   ./configure -openmm -openblas -fftw3 -zlib -bzlib -shared -openmp -arpack gnu
fi

make install

rsync -a README.md config.h LICENSE dat bin lib $PREFIX

# mkdir -p $PREFIX/doc
# rsync -a doc/cpptraj.pdf $PREFIX/doc

# Export CPPTRAJHOME automatically
mkdir -p ${PREFIX}/etc/conda/{activate,deactivate}.d
cp ${RECIPE_DIR}/activate.sh ${PREFIX}/etc/conda/activate.d/cpptraj.sh
cp ${RECIPE_DIR}/deactivate.sh ${PREFIX}/etc/conda/deactivate.d/cpptraj.sh
