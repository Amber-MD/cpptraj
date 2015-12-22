#!/bin/bash

. ../MasterTest.sh

CleanFiles rotate.in fromMatrices.crd

INPUT="-i rotate.in"
cat > rotate.in <<EOF
parm ../tz2.truncoct.parm7
parmstrip :WAT
trajin ../Test_RMSD/tz2.norotate.crd.save
readdata ../Test_RMSD/rmatrices.dat.save name RM mat3x3
rotate usedata RM
# Note: Not doing any translations
#rms reftraj tz2.rotate.crd.save out fromMatrices.dat
trajout fromMatrices.crd
EOF
RunCpptraj "Rotation of coords from matrices"
DoTest fromMatrices.crd.save fromMatrices.crd

EndTest
exit 0
