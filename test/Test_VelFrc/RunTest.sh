#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in CrdFrcVel.nc Vel.crd Frc.crd Vel1.crd Frc1.crd \
           CrdFrcVel.ncrst.? fromncrst.nc

TESTNAME='Read separate velocity/force trajectory data tests'
Requires notparallel netcdf

INPUT="-i cpptraj.in"

cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
trajin short.crd mdvel short.vel mdfrc short.frc
trajout CrdFrcVel.nc
EOF
RunCpptraj "Test combining coordinates, velocities, and forces."

cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
trajin CrdFrcVel.nc usevelascoords
trajout Vel.crd
EOF
RunCpptraj "Test using velocities as coordinates."
DoTest Vel.crd.save Vel.crd

cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
trajin CrdFrcVel.nc
trajout Vel1.crd mdvel
EOF
RunCpptraj "Test writing velocities (MDVEL)"
DoTest Vel.crd.save Vel1.crd

cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
trajin CrdFrcVel.nc usefrcascoords
trajout Frc.crd
EOF
RunCpptraj "Test using forces as coordinates."
DoTest Frc.crd.save Frc.crd

cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
trajin CrdFrcVel.nc
trajout Frc1.crd mdfrc
EOF
RunCpptraj "Test writing forces (MDFRC)"
DoTest Frc.crd.save Frc1.crd

cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
trajin short.crd mdvel short.vel mdfrc short.frc
trajout CrdFrcVel.ncrst
EOF
RunCpptraj "Test writing combined coords/velocity/force NetCDF restart."

cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
trajin CrdFrcVel.ncrst.1
trajin CrdFrcVel.ncrst.2
trajout fromncrst.nc
EOF
RunCpptraj "Test reading combined coords/velocity/force NetCDF restart."
NcTest CrdFrcVel.nc fromncrst.nc

EndTest
exit 0
