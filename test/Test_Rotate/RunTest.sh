#!/bin/bash

. ../MasterTest.sh

CleanFiles rotate.in fromMatrices.crd TCS.rotated.mol2 inverse.crd \
           tz2.rotate.rst7 tz2.*.rotate.rst7 tz2.?.rst7 matrices.dat \
           rotations.dat

TESTNAME='Coordinates rotation tests'
Requires maxthreads 10

INPUT="-i rotate.in"
cat > rotate.in <<EOF
parm ../tz2.truncoct.parm7
parmstrip :WAT
trajin ../Test_RMSD/tz2.norotate.crd.save
readdata ../Test_RMSD/rmatrices.dat.save name RM mat3x3
rotate usedata RM
# Note: Not doing any translations
#rms reftraj tz2.rotate.crd.save out fromMatrices.dat
outtraj fromMatrices.crd
rotate usedata RM inverse
outtraj inverse.crd
EOF
RunCpptraj "Rotation (and inverse) of coords from matrices"
DoTest fromMatrices.crd.save fromMatrices.crd
DoTest ../Test_RMSD/tz2.norotate.crd.save inverse.crd

UNITNAME='Rotate with defined axis'
CheckFor netcdf maxthreads 1
if [ $? -eq 0 ] ; then
  cat > rotate.in <<EOF
parm ../FtuFabI.NAD.TCL.parm7
trajin ../FtuFabI.NAD.TCL.nc 1 1
rotate :270 axis0 :270@C1,C2,C3,C4,C5,C6 axis1 :270@C7,C8,C9,C10,C11,C12 90.0
strip :1-269
trajout TCS.rotated.mol2
EOF
  RunCpptraj "$UNITNAME"
  DoTest TCS.rotated.mol2.save TCS.rotated.mol2
fi

UNITNAME='Basic rotation around axes'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > rotate.in <<EOF
parm ../tz2.parm7

trajin ../tz2.rst7
center origin
#outtraj tz2.0.rst7
rotate * x 30
#outtraj tz2.1.rst7
rotate * y 45
#outtraj tz2.2.rst7
rotate * z 60
trajout tz2.separate.rotate.rst7
run

clear trajin

trajin ../tz2.rst7
center origin
#debug actions 1
rotate * x 30 y 45 z 60
trajout tz2.combined.rotate.rst7
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2.separate.rotate.rst7.save tz2.separate.rotate.rst7
  DoTest tz2.separate.rotate.rst7 tz2.combined.rotate.rst7
  #DoTest tz2.rotate.rst7.save tz2.rotate.rst7
fi

UNITNAME='Calculate rotations from rotation matrices test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > rotate.in <<EOF
parm ../tz2.parm7
# Want to go from the original coords (target) to rotated coords (reference)
reference tz2.separate.rotate.rst7.save name REF
trajin ../tz2.rst7
rms R0 reference savematrices matricesout matrices.dat
rotate calcfrom R0[RM] name Rot out rotations.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest rotations.dat.save rotations.dat
fi

EndTest
exit 0
