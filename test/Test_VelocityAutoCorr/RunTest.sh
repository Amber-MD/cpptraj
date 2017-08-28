#!/bin/bash

. ../MasterTest.sh

CleanFiles vac.in VAC.agr VAC2.dat diff.dat

INPUT='-i vac.in'

cat > vac.in <<EOF
parm ../Test_systemVF/systemVF.parm7
trajin ../Test_systemVF/systemVF.nc
velocityautocorr out VAC.agr  Vel           norm :WAT@O diffout diff.dat
velocityautocorr out VAC.agr  Crd           norm :WAT@O usecoords
velocityautocorr out VAC2.dat Direct direct norm :WAT@O
EOF
RunCpptraj "Velocity autocorrelation test"
DoTest VAC.agr.save VAC.agr
DoTest diff.dat.save diff.dat
DoTest VAC2.dat.save VAC2.dat

EndTest
exit 0
