#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in surf.dat tsurf.dat ral.surf.dat

INPUT="-i cpptraj.in"

# Test 1
UNITNAME='Surface calculation test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc
surf All out surf.dat
EOF
  TOP="../DPDP.parm7"
  RunCpptraj "$UNITNAME"
  DoTest surf.dat.save surf.dat
fi

# Test 2
UNITNAME='Partial Surface calculation test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc
surf R1-12 out tsurf.dat :1-12
EOF
  TOP="../DPDP.parm7"
  RunCpptraj "$UNITNAME"
  DoTest tsurf.dat.save tsurf.dat
fi

# Test 3
UNITNAME='LCPO test with GAFF atom types'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  TOP=""
  cat > cpptraj.in <<EOF
parm RAL.sol.top
trajin RAL.crd
surf S0 out ral.surf.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest ral.surf.dat.save ral.surf.dat
fi

EndTest

exit 0
