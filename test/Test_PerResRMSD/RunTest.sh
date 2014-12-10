#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles PerResRMSD.agr cpptraj.in test.dat perresavg.dat center.agr

# Test 1
CheckNetcdf
cat > cpptraj.in <<EOF
parm ../tz2.truncoct.parm7
reference ../tz2.truncoct.nc 1
trajin ../tz2.truncoct.nc  
rmsd :2-11 refindex 0 perres perresout PerResRMSD.agr
EOF
INPUT="-i cpptraj.in"
RunCpptraj "Per-Residue RMSD Test."
DoTest PerResRMSD.agr.save PerResRMSD.agr
CheckTest

# Test 2
cat > cpptraj.in <<EOF
parm ../tz2.truncoct.parm7
reference ../tz2.truncoct.nc 1
trajin ../tz2.truncoct.nc 2 4 
rmsd :2-11 refindex 0 perres perresavg perresavg.dat
EOF
INPUT="-i cpptraj.in"
RunCpptraj "Per-Residue RMSD Test with averaging."
DoTest perresavg.dat.save perresavg.dat
CheckTest

# Test 3
cat > cpptraj.in <<EOF
parm ../tz2.truncoct.parm7
reference ../tz2.truncoct.nc 1
trajin ../tz2.truncoct.nc  
rmsd :2-11 refindex 0 perres perresout center.agr range 1 perrescenter 
EOF
INPUT="-i cpptraj.in"
RunCpptraj "Per-Residue RMSD Test with residue centering."
DoTest center.agr.save center.agr
CheckTest

EndTest

exit 0
