#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles radgyr.in radgyr.dat radgyr.mass.dat 

# Test 1
CheckNetcdf
cat > radgyr.in <<EOF
noprogress
parm ../ChainA-tip3p.parm7
trajin ../run0.nc 1 10
radgyr out radgyr.dat :1-268
radgyr out radgyr.mass.dat :1-268 mass
EOF
INPUT="-i radgyr.in"
RunCpptraj "Radius of gyration command test."
DoTest radgyr.dat.save radgyr.dat
DoTest radgyr.mass.dat.save radgyr.mass.dat 

CheckTest

EndTest

exit 0
