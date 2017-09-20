#!/bin/bash

. ../MasterTest.sh

CleanFiles min.in min.dat

TESTNAME='Minimum non-self imaged distance test'
Requires netcdf maxthreads 10

INPUT="-i min.in"

cat > min.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
minimage m1 :1-13 :1-13 out min.dat
minimage m2 :1-13 :1-13 maskcenter out min.dat
EOF
RunCpptraj "Minimum Image test."
DoTest min.dat.save min.dat
EndTest
exit 0
