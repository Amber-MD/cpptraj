#!/bin/bash

. ../MasterTest.sh

CleanFiles ws.in ws.agr

CheckNetcdf
TOP=../tz2.truncoct.parm7
INPUT="ws.in"
cat > ws.in <<EOF
trajin ../tz2.truncoct.nc
watershell !:WAT ws.agr
EOF
RunCpptraj "Watershell Test"
DoTest ws.agr.save ws.agr
EndTest

exit 0
