#!/bin/bash

. ../MasterTest.sh

CleanFiles cc.in cc.gnu dist.agr
RequiresNetcdf "CrossCorr test"
INPUT="-i cc.in"
cat > cc.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
distance d1-13 :1 :13
distance d3-11 :3 :11
distance d5-9  :5 :9
crosscorr name CC d1-13 d3-11 d5-9 out cc.gnu
EOF
RunCpptraj "CrossCorr test."
DoTest cc.gnu.save cc.gnu
EndTest

exit 0
