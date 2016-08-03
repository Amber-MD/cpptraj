#!/bin/bash

. ../MasterTest.sh

CleanFiles vac.in VAC.agr

INPUT='-i vac.in'

cat > vac.in <<EOF
parm ../Test_systemVF/systemVF.parm7
trajin ../Test_systemVF/systemVF.nc
velocityautocorr out VAC.agr Vel usevelocity norm :WAT@O
velocityautocorr out VAC.agr Crd norm :WAT@O
EOF
RunCpptraj "Velocity autocorrelation test"
DoTest VAC.agr.save VAC.agr

EndTest
exit 0
