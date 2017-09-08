#!/bin/bash

. ../MasterTest.sh

CleanFiles vol.in vol.dat

TESTNAME='Volume tests'
Requires netcdf

INPUT='-i vol.in'

cat > vol.in <<EOF
parm ../tz2.truncoct.parm7 [OCT]
parm ../tz2.ortho.parm7 [ORTHO]

loadcrd ../tz2.truncoct.nc name T1 parm [OCT]
loadcrd ../tz2.ortho.nc name T2 parm [ORTHO]

crdaction T1 volume Oct out vol.dat
crdaction T2 volume Ortho out vol.dat
EOF
RunCpptraj "$TESTNAME"
DoTest vol.dat.save vol.dat

EndTest
exit 0
