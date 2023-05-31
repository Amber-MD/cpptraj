#!/bin/bash

. ../MasterTest.sh

CleanFiles tordiff.in frc.dat tor.dat

TESTNAME='Toroidal-view-preserving diffusion calculation tests'
Requires netcdf maxthreads 1 

INPUT='-i tordiff.in'

UNITNAME='Basic toroidal-view-preserving diffusion calculation test'
cat > tordiff.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc

diffusion FRC :WAT@O out frc.dat
tordiff TOR :WAT@O out tor.dat
EOF
RunCpptraj "$UNITNAME"
DoTest tor.dat.save tor.dat

EndTest
exit 0
