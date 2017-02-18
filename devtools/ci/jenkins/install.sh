# This script is sourced from the main build script, and the executing directory
# is the top-level cpptraj directory. The environment variable COMPILER_FLAGS is
# set, as is label. This is set up exclusively to test the Intel compilers.

# Load the Intel compilers (this also sets MKL_HOME)
export TEST_OS="${label}"
if [ "${label}" = "linux" ]; then
  module load intel openmpi-intel amber/17 cuda
  # Export the shader model of a GTX-680 (which is on the Jenkins machine)
  export SHADER_MODEL=sm_30
  ./configure ${COMPILER_FLAGS} -mkl intel
else
  # Mac OS X
  ./configure -macAccelerate --with-fftw3=/opt/local --with-netcdf=/opt/local -noarpack ${COMPILER_FLAGS} clang
fi

# If the compiler flag is -mpi, set DO_PARALLEL. Test both 2 and 4 CPUs for
# parallel builds.
make -j4 install
if [ "${COMPILER_FLAGS}" = "-mpi" ]; then
  export DO_PARALLEL='mpirun -np 2'
  make check
  export DO_PARALLEL='mpirun -np 4'
  make check
elif [ "${COMPILER_FLAGS}" = "-openmp" ]; then
  export OPT=openmp OMP_NUM_THREADS=2
  make check
  export OMP_NUM_THREADS=4
  make check
elif [ "${COMPILER_FLAGS}" = "-cuda" ]; then
  export OPT=cuda
  make check
else
  make check
fi
