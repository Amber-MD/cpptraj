#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd.dat rms.mass.in rmsd.mass.dat rms.reftraj.in rmsd.reftraj.dat

CheckNetcdf
TOP="../tz2.truncoct.parm7"
INPUT="rms.in"

# Test rmsd, mass-weighted rmsd, rmsd to reference traj.
cat > rms.in <<EOF
noprogress
trajin ../tz2.truncoct.nc

rms Res2-11 first :2-11 out rmsd.dat
rms Res2-11_mass first :2-11 out rmsd.mass.dat mass
rms Res2_11_traj reftraj ../tz2.truncoct.nc :2-11 out rmsd.reftraj.dat

EOF
RunCpptraj "RMSD Tests."
DoTest rmsd.dat.save rmsd.dat
DoTest rmsd.mass.dat.save rmsd.mass.dat
DoTest rmsd.reftraj.dat.save rmsd.reftraj.dat

EndTest

exit 0
