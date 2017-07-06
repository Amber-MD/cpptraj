#!/bin/bash

. ../MasterTest.sh

CleanFiles vol.in vol.dat

CheckNetcdf

INPUT='-i vol.in'

cat > vol.in <<EOF
parm ../tz2.truncoct.parm7 [OCT]
parm ../tz2.ortho.parm7 [ORTHO]

loadcrd ../tz2.truncoct.nc name T1 parm [OCT]
loadcrd ../tz2.ortho.nc name T2 parm [ORTHO]

crdaction T1 volume Oct out vol.dat
crdaction T2 volume Ortho out vol.dat
crdaction T1 volume density dOct out vol.dat
crdaction T2 volume density dOrtho out vol.dat
EOF
RunCpptraj "Volume/density tests."
DoTest vol.dat.save vol.dat

EndTest
exit 0
