#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles jcoupling.in Jcoupling.dat

# Test 1
CheckNetcdf
cat > jcoupling.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc 1 1
jcoupling outfile Jcoupling.dat
EOF
INPUT="-i jcoupling.in"
RunCpptraj "J-Coupling command test."
DoTest Jcoupling.dat.save Jcoupling.dat

CheckTest

EndTest

exit 0
