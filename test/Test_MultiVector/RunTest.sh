#!/bin/bash

. ../MasterTest.sh

CleanFiles mv.in NH.dat V2.dat
CheckNetcdf
INPUT="-i mv.in"
cat > mv.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
multivector NH name1 N name2 H ired out NH.dat
vector V2 :2@N :2@H out V2.dat
EOF
RunCpptraj "Multivector test."
DoTest NH.dat.save NH.dat
EndTest
exit 0
