#!/bin/bash

. ../MasterTest.sh
CleanFiles corr.in corr.dat cross.dat
TESTNAME='Correlation test'
Requires netcdf
INPUT="-i corr.in"
cat > corr.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc

distance d1 :2 :12
radgyr rg :1-3 mass nomax 

corr d1 d1 out corr.dat 
corr d1 rg out cross.dat 
EOF
RunCpptraj "$TESTNAME"
DoTest corr.dat.save corr.dat
DoTest cross.dat.save cross.dat
EndTest

exit 0
