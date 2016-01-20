#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in surf.dat tsurf.dat ral.surf.dat

CheckNetcdf
INPUT="-i cpptraj.in"
# Test 1
cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc
surf All out surf.dat
EOF
TOP="../DPDP.parm7"
RunCpptraj "Surface calculation test."
DoTest surf.dat.save surf.dat

# Test 2
cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc
surf R1-12 out tsurf.dat :1-12
EOF
TOP="../DPDP.parm7"
RunCpptraj "Partial Surface calculation test."
DoTest tsurf.dat.save tsurf.dat

# Test 3
MaxThreads 1 "LCPO test with GAFF atom types."
if [[ $? -eq 0 ]] ; then
  TOP=""
  cat > cpptraj.in <<EOF
parm RAL.sol.top
trajin RAL.crd
surf S0 out ral.surf.dat
EOF
  RunCpptraj "LCPO test with GAFF atom types."
  DoTest ral.surf.dat.save ral.surf.dat
fi

EndTest

exit 0
