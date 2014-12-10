#!/bin/bash

. ../MasterTest.sh

CleanFiles rotdif.in rvecs.dat matrices.dat deffs.dat rotdif.out
CheckNetcdf
CheckPtrajAnalyze

INPUT="-i rotdif.in"
cat > rotdif.in <<EOF
parm ../tz2.parm7
reference avgstruct.pdb [tz2avg] 
trajin ../tz2.nc 
rotdif rseed 1 nvecs 10 ref [tz2avg] @CA,C,N,O dt 0.002 tf 0.190 \
       itmax 500 tol 0.000001 d0 0.03 order 2 rvecout rvecs.dat \
       rmout matrices.dat deffout deffs.dat outfile rotdif.out
EOF
RunCpptraj "Rotdif test"
DoTest rvecs.dat.save rvecs.dat
DoTest matrices.dat.save matrices.dat
DoTest deffs.dat.save deffs.dat
DoTest rotdif.out.save rotdif.out

CheckTest
EndTest

exit 0

