#!/bin/bash

. ../MasterTest.sh

CleanFiles ac.in ac.agr dist.agr vac.agr
TESTNAME='AutoCorr test'
RequiresNetcdf "$TESTNAME"
INPUT="-i ac.in"
cat > ac.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
distance d1-13 :1 :13
distance d3-11 :3 :11
distance d5-9  :5 :9
vector v1 :2@N :2@H
#create dist.agr d1-13 d3-11 d5-9
autocorr name AC d1-13 d3-11 d5-9 v1 out ac.agr lagmax 50
EOF
RunCpptraj "AutoCorr test."
DoTest ac.agr.save ac.agr
EndTest

exit 0
