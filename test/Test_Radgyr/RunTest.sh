#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles radgyr.in radgyr.dat radgyr.mass.dat 

# Test 1
CheckNetcdf
cat > radgyr.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 10
radgyr Res1-13 out radgyr.dat :1-13
radgyr out radgyr.mass.dat :1-13 mass
EOF
INPUT="-i radgyr.in"
RunCpptraj "Radius of gyration command test."
DoTest radgyr.dat.save radgyr.dat
DoTest radgyr.mass.dat.save radgyr.mass.dat 

CheckTest

EndTest

exit 0
