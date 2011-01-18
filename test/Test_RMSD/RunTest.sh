#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd.dat rms.mass.in rmsd.mass.dat

CheckNetcdf
TOP="../ChainA-tip3p.parm7"

# Test 1
cat > rms.in <<EOF
trajin ../run0.nc
rms first :10-260 out rmsd.dat
EOF
INPUT="rms.in"
RunCpptraj "RMSD Test."
DoTest rmsd.dat.save rmsd.dat

# Test 2
cat > rms.mass.in <<EOF
trajin ../run0.nc
rms first :10-260 out rmsd.mass.dat mass
EOF
INPUT="rms.mass.in"
RunCpptraj "RMSD test with mass weighting."
DoTest rmsd.mass.dat.save rmsd.mass.dat
CheckTest

EndTest

exit 0
