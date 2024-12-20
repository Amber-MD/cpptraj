#!/bin/sh

./configure -openblas -fftw3 -zlib -bzlib -shared gnu

make install

rsync -a README.md Makefile config.h LICENSE dat bin lib $PREFIX

# mkdir -p $PREFIX/doc
# rsync -a doc/cpptraj.pdf $PREFIX/doc

# Export CPPTRAJHOME automatically
mkdir -p ${PREFIX}/etc/conda/{activate,deactivate}.d
cp ${RECIPE_DIR}/activate.sh ${PREFIX}/etc/conda/activate.d/cpptraj.sh
cp ${RECIPE_DIR}/deactivate.sh ${PREFIX}/etc/conda/deactivate.d/cpptraj.sh
