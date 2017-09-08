#!/bin/bash

. ../MasterTest.sh

CleanFiles molsurf.in msurf.dat
TESTNAME='Molsurf test'
Requires netcdf
# Molsurf test
cat > molsurf.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
molsurf GB    out msurf.dat
molsurf PARSE out msurf.dat radii parse
molsurf VDW   out msurf.dat radii vdw
EOF
INPUT="-i molsurf.in"
RunCpptraj "$TESTNAME"
DoTest msurf.dat.save msurf.dat
EndTest

exit 0
