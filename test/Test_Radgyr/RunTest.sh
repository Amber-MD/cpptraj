#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles radgyr.in radgyr.dat radgyr.mass.dat 

TESTNAME='Radius of gyration command test'
Requires netcdf maxthreads 10

# Test 1
cat > radgyr.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 10
radgyr Res1-13 out radgyr.dat :1-13
radgyr Res1-13_mass out radgyr.mass.dat :1-13 mass
EOF
INPUT="-i radgyr.in"
RunCpptraj "$TESTNAME"
DoTest radgyr.dat.save radgyr.dat
DoTest radgyr.mass.dat.save radgyr.mass.dat 

EndTest

exit 0
