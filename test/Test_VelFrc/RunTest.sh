#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in CrdFrcVel.nc Vel.crd Frc.crd

NotParallel "Separate velocity/force"
if [ "$?" -eq 1 ] ; then
  EndTest
  exit 0
fi
CheckNetcdf

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
trajin CrdFrcVel.nc usefrcascoords
trajout Frc.crd
EOF
RunCpptraj "Test using forces as coordinates."
DoTest Frc.crd.save Frc.crd

EndTest
exit 0
