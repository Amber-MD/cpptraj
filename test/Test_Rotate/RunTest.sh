#!/bin/bash

. ../MasterTest.sh

CleanFiles rotate.in fromMatrices.crd TCS.rotated.mol2

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

MaxThreads 1 "Rotate with defined axis."
if [[ $? -eq 0 ]] ; then
  cat > rotate.in <<EOF
parm ../FtuFabI.NAD.TCL.parm7
trajin ../FtuFabI.NAD.TCL.nc 1 1
rotate :270 axis0 :270@C1,C2,C3,C4,C5,C6 axis1 :270@C7,C8,C9,C10,C11,C12 90.0
strip :1-269
trajout TCS.rotated.mol2
EOF
  RunCpptraj "Rotate with defined axis."
  DoTest TCS.rotated.mol2.save TCS.rotated.mol2
fi
EndTest
exit 0
