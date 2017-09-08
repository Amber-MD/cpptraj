#!/bin/bash

. ../MasterTest.sh

CleanFiles rotdif.in rvecs.dat matrices.dat deffs.dat rotdif.out
TESTNAME='Rotational diffusion calculation test'
Requires netcdf mathlib

INPUT="-i rotdif.in"
cat > rotdif.in <<EOF
parm ../tz2.parm7
reference avgstruct.pdb [tz2avg] 
trajin ../tz2.nc
rms R0 ref [tz2avg] @CA,C,N,O savematrices 
rotdif rmatrix R0[RM] rseed 1 nvecs 10 dt 0.002 tf 0.190 ncorr 101 \
       itmax 500 tol 0.000001 d0 0.03 order 2 rvecout rvecs.dat \
       rmout matrices.dat deffout deffs.dat outfile rotdif.out
EOF
RunCpptraj "$TESTNAME"
DoTest rvecs.dat.save rvecs.dat
DoTest matrices.dat.save matrices.dat
DoTest deffs.dat.save deffs.dat
DoTest rotdif.out.save rotdif.out

EndTest

exit 0

