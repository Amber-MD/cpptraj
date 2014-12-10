#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in surf.dat tsurf.dat

# Test 1
CheckNetcdf
cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc
surf out surf.dat 
EOF
INPUT="-i cpptraj.in"
TOP="../DPDP.parm7"
RunCpptraj "Surface calculation test."
DoTest surf.dat.save surf.dat
CheckTest

# Test 2
cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc
surf out tsurf.dat :1-12 
EOF
INPUT="-i cpptraj.in"
TOP="../DPDP.parm7"
RunCpptraj "Partial Surface calculation test."
DoTest tsurf.dat.save tsurf.dat
CheckTest

EndTest

exit 0
