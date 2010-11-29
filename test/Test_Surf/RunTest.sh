#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in surf.dat

# Test 1
CheckNetcdf
cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc
#trajout dpdp.short.nc netcdf
surf out surf.dat 
EOF
INPUT="-i cpptraj.in"
TOP="../DPDP.mod.GA12.parm7"
RunCpptraj "Surface calculation test."
DoTest surf.dat.save surf.dat
CheckTest

EndTest

exit 0
