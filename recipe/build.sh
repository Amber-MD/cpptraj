#!/bin/sh

./configure -openmp -openblas -shared gnu
make install

rsync -a README.md LICENSE bin lib $PREFIX

# mkdir -p $PREFIX/doc
# rsync -a doc/cpptraj.pdf $PREFIX/doc

# Export AMBERCLASSICHOME automatically
# mkdir -p ${PREFIX}/etc/conda/{activate,deactivate}.d
# cp ${RECIPE_DIR}/activate.sh ${PREFIX}/etc/conda/activate.d/amberclassic.sh
# cp ${RECIPE_DIR}/deactivate.sh ${PREFIX}/etc/conda/deactivate.d/amberclassic.sh
