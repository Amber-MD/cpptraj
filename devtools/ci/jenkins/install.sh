# This script is sourced from the main build script, and the executing directory
# is the top-level cpptraj directory. The environment variable COMPILER_FLAGS is
# set. This is set up exclusively to test the Intel compilers.

# Load the Intel compilers (this also sets MKL_HOME)
module load intel openmpi-intel

./configure ${COMPILER_FLAGS} -mkl intel

# If the compiler flag is -mpi, set DO_PARALLEL. Test both 2 and 4 CPUs for
# parallel builds.
make -j6 install
if [ "${COMPILER_FLAGS}" = "-mpi" ]; then
  export DO_PARALLEL='mpirun -np 2'
  make check
  export DO_PARALLEL='mpirun -np 4'
  make check
elif [ "${COMPILER_FLAGS}" = "-openmp" ]; then
  export OMP_NUM_THREADS=2
  make check
  export OMP_NUM_THREADS=4
  make check
else
  make check
fi
