#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in surf.dat tsurf.dat

# Test 1
CheckNetcdf
cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc
surf All out surf.dat
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
surf R1-12 out tsurf.dat :1-12
EOF
INPUT="-i cpptraj.in"
TOP="../DPDP.parm7"
RunCpptraj "Partial Surface calculation test."
DoTest tsurf.dat.save tsurf.dat
CheckTest

# Test 3
INPUT="-i cpptraj.in"
TOP=""
cat > cpptraj.in <<EOF
parm RAL.sol.top
trajin RAL.crd
surf S0 out ral.surf.dat
EOF
RunCpptraj "LCPO test with GAFF atom types."
DoTest ral.surf.dat.save ral.surf.dat


EndTest

exit 0
