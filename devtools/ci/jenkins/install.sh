# MKL_HOME gets set by the Jenkinsfile even if it's not actually installed. So unset
# it if it's set and not actually a directory
test -d $MKL_HOME || unset MKL_HOME

# This script is sourced from the main build script, and the executing directory
# is the top-level cpptraj directory. The environment variable COMPILER_FLAGS is
# set, as is label. This is set up exclusively to test the Intel compilers.

# Debug -- print environment
env

if [ "${OPERATING_SYSTEM}" = "linux" ]; then
  bash +e -x ./configure --with-fftw3 --with-netcdf ${COMPILER_FLAGS} ${COMPILER}
else
  # Mac OS X
  bash +e ./configure -macAccelerate --with-fftw3=/opt/local --with-netcdf=/opt/local -noarpack ${COMPILER_FLAGS} clang
fi

compiler_flags_contains() {
  if [ "x`echo "${COMPILER_FLAGS}" | sed -e "s/$1//g"`" = "x${COMPILER_FLAGS}" ]; then
    echo "no" # It's the same, so $1 cannot be in COMPILER_FLAGS
  else
    echo "yes"
  fi
}

# If the compiler flag is -mpi, set DO_PARALLEL. Test both 2 and 4 CPUs for
# parallel builds.
make -j4 install

echo "cpptraj build is complete!"

if [ `compiler_flags_contains -mpi` = "yes" ]; then
  export DO_PARALLEL='mpirun -np 2'
  make check
  export DO_PARALLEL='mpirun -np 4'
  make check
elif [ `compiler_flags_contains -openmp` = "yes" ]; then
  export OPT=openmp OMP_NUM_THREADS=2
  make check
  export OMP_NUM_THREADS=4
  make check
elif [ `compiler_flags_contains -cuda` = "yes" ]; then
  export OPT=cuda
  make check
else
  make check
fi
