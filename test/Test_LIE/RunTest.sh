#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles lie.in TCL.out

TESTNAME='LIE test'
Requires netcdf maxthreads 10

INPUT="lie.in"
TOP=../FtuFabI.NAD.TCL.parm7
cat > lie.in <<EOF
trajin ../FtuFabI.NAD.TCL.nc
lie LIE :TCS out TCL.out cutvdw 12 cutelec 12
EOF
RunCpptraj "$TESTNAME"
DoTest TCL.out.save TCL.out

EndTest

exit 0
