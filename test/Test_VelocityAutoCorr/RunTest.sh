#!/bin/bash

. ../MasterTest.sh

CleanFiles vac.in VAC.agr VAC2.dat diff.dat

TESTNAME='Velocity autocorrelation test'
Requires netcdf maxthreads 5

INPUT='-i vac.in'

cat > vac.in <<EOF
parm ../Test_systemVF/systemVF.parm7
trajin ../Test_systemVF/systemVF.nc
velocityautocorr out VAC.agr Vel usevelocity norm :WAT@O diffout diff.dat
velocityautocorr out VAC.agr Crd norm :WAT@O
velocityautocorr out VAC2.dat Direct usevelocity direct norm :WAT@O
EOF
RunCpptraj "$TESTNAME"
DoTest VAC.agr.save VAC.agr
DoTest diff.dat.save diff.dat
DoTest VAC2.dat.save VAC2.dat

EndTest
exit 0
