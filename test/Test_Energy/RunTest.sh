#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in ene.agr short.dat

INPUT="-i ene.in"

TESTNAME='Simple energy tests'
Requires netcdf

cat > ene.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
energy DPDP out ene.agr
EOF
RunCpptraj "$TESTNAME"
DoTest ene.agr.save ene.agr

UNITNAME='Test kinetic energy calculation'
CheckFor maxthreads 2
if [ $? -eq 0 ] ; then
  cat > ene.in <<EOF
parm ../tz2.nhe.parm7
trajin ../Test_VelFrc/short.crd mdvel ../Test_VelFrc/short.vel mdfrc ../Test_VelFrc/short.frc
energy kinetic Short out short.dat dt 0.002
EOF
  RunCpptraj "$UNITNAME"
  DoTest short.dat.save short.dat
fi
EndTest
exit 0
