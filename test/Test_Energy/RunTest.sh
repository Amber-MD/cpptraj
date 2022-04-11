#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in ene.agr short.dat tz2.dat strip.dat elec.dat vdw.dat

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

UNITNAME='Test total energy calculation'
CheckFor maxthreads 2
if [ $? -eq 0 ] ; then
  cat > ene.in <<EOF
parm ../tz2.nhe.parm7
trajin ../Test_VelFrc/short.crd mdvel ../Test_VelFrc/short.vel mdfrc ../Test_VelFrc/short.frc
energy Tz2 out tz2.dat dt 0.002
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2.dat.save tz2.dat
fi

UNITNAME='Test partial energy calculation'
CheckFor maxthreads 10
if [ $? -eq 0 ] ; then
  cat > ene.in <<EOF
parm ../FtuFabI.NAD.TCL.parm7
trajin ../FtuFabI.NAD.TCL.nc
energy out strip.dat :NDP ENE1
energy out elec.dat  :NDP ENE2 elec
energy out vdw.dat   :NDP ENE3 vdw
EOF
  RunCpptraj "$UNITNAME"
  DoTest strip.dat.save strip.dat
  DoTest elec.dat.save elec.dat
  DoTest vdw.dat.save vdw.dat
fi

EndTest
exit 0
