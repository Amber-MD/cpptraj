#!/bin/bash

. ../MasterTest.sh

CleanFiles diff.in

INPUT='diff.in'

TESTNAME='Diffusion analysis tests'
Requires netcdf notparallel

cat > diff.in <<EOF
parm ../tz2.ortho.parm7
loadcrd ../tz2.ortho.nc name TZ2 1 10 
crdaction TZ2 unwrap bymol
runanalysis calcdiffusion crdset TZ2 out tz2.diff.dat :WAT@O
EOF
RunCpptraj "Diffusion analysis calculation test"
DoTest tz2.diff.dat.save tz2.diff.dat

EndTest
exit 0
