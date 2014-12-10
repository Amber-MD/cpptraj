#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd.dat rms.mass.in rmsd.mass.dat rms.reftraj.in rmsd.reftraj.dat

CheckNetcdf
TOP="../tz2.truncoct.parm7"

# Test 1
cat > rms.in <<EOF
noprogress
trajin ../tz2.truncoct.nc
rms first :2-11 out rmsd.dat
EOF
INPUT="rms.in"
RunCpptraj "RMSD Test."
DoTest rmsd.dat.save rmsd.dat

# Test 2
cat > rms.mass.in <<EOF
noprogress
trajin ../tz2.truncoct.nc
rms first :2-11 out rmsd.mass.dat mass
EOF
INPUT="rms.mass.in"
RunCpptraj "RMSD test with mass weighting."
DoTest rmsd.mass.dat.save rmsd.mass.dat

# Test 3
cat > rms.reftraj.in <<EOF
noprogress
trajin ../tz2.truncoct.nc
rms reftraj ../tz2.truncoct.nc :2-11 out rmsd.reftraj.dat
EOF
INPUT="rms.reftraj.in"
RunCpptraj "RMSD test with reference trajectory."
DoTest rmsd.reftraj.dat.save rmsd.reftraj.dat

CheckTest

EndTest

exit 0
