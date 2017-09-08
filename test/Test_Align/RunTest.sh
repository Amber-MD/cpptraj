#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd.dat rmsd.mass.dat rmsd.reftraj.dat

TESTNAME='Align tests'
Requires netcdf
TOP="../tz2.truncoct.parm7"
INPUT="rms.in"

# Test rmsd, mass-weighted rmsd, rmsd to reference traj.
cat > rms.in <<EOF
noprogress
trajin ../tz2.truncoct.nc
loadtraj ../tz2.truncoct.nc name TZ2

align :2-11 first
rms Res2-11 nofit :2-11 out rmsd.dat
align :2-11 mass first
rms Res2-11_mass nofit :2-11 out rmsd.mass.dat mass
align :2-11 reftraj TZ2
rms Res2_11_traj nofit reftraj TZ2 :2-11 out rmsd.reftraj.dat
EOF
RunCpptraj "Align Tests."
DoTest ../Test_RMSD/rmsd.dat.save rmsd.dat
DoTest ../Test_RMSD/rmsd.mass.dat.save rmsd.mass.dat
DoTest ../Test_RMSD/rmsd.reftraj.dat.save rmsd.reftraj.dat

EndTest

exit 0
