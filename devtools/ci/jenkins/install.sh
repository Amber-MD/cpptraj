# This script is sourced from the main build script, and the executing directory
# is the top-level cpptraj directory. The environment variable COMPILER_FLAGS is
# set. This is set up exclusively to test the Intel compilers.

# Load the Intel compilers (this also sets MKL_HOME). It also assumes that if
# -mpi is set, so is DO_PARALLEL
module load intel openmpi-intel

./configure ${COMPILER_FLAGS} -mkl intel
make -j6 install
make check
