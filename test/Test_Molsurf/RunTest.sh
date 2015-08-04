#!/bin/bash

. ../MasterTest.sh

CleanFiles molsurf.in msurf.dat

# Molsurf test
cat > molsurf.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
molsurf out msurf.dat
#trajout pqr pdb dumpq multi
EOF
INPUT="-i molsurf.in"
RunCpptraj "Molsurf test."
DoTest msurf.dat.save msurf.dat
CheckTest
EndTest

exit 0
