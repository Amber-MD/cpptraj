#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd.dat rmsd.mass.dat rmsd.reftraj.dat report.?.dat

TESTNAME='Align tests'
Requires netcdf maxthreads 10
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

# Check before and after alignment should be the same
cat > rms.in <<EOF
trajin ../tz2.truncoct.nc
check reportfile report.0.dat
center :1-13
align first :1-13@CA
image triclinic :WAT
check reportfile report.1.dat
EOF
RunCpptraj "Align, unit cell rotation/imaging test."
DoTest report.0.dat report.1.dat

EndTest

exit 0
