#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles lie.in TCL.out

INPUT="lie.in"
TOP=../FtuFabI.NAD.TCL.parm7
# Check libraries
CheckNetcdf
cat > lie.in <<EOF
trajin ../FtuFabI.NAD.TCL.nc
lie LIE :TCS out TCL.out cutvdw 12 cutelec 12
EOF
RunCpptraj "LIE test, TCL"
DoTest TCL.out.save TCL.out

EndTest

exit 0
