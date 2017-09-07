#!/bin/bash

. ../MasterTest.sh

CleanFiles runavg.in running_avg.dat cumulative_avg.dat distances.dat
TESTNAME='Analysis Running Average'
RequiresNetcdf "$TESTNAME"
INPUT="runavg.in"
TOP="../tz2.parm7"
cat > $INPUT <<EOF
trajin ../tz2.nc 1 14
distance d1 :1 :7 out distances.dat
distance d2 :10 :13 out distances.dat
analyze runningavg d1 d2 out cumulative_avg.dat cumulative
analyze runningavg d1 d2 out running_avg.dat window 10
EOF
RunCpptraj "$TESTNAME"
DoTest running_avg.dat.save running_avg.dat
DoTest cumulative_avg.dat.save cumulative_avg.dat
EndTest

exit 0
