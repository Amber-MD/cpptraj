#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in ene.dat ene1.dat long.dat strip.dat directsum.0 ewald.dat \
           ew_tz2.dat ew_tz2_10.dat tz2_ortho.dat directsum.0 run9.dat pme.nacl.dat \
            ew_tz2o.dat

INPUT="-i ene.in"
TESTNAME='Particle mesh Ewald tests'
Requires libpme maxthreads 10

Simple() {
  UNITNAME='Particle mesh Ewald test (simple)'
  CheckFor maxthreads 1
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm test.mol2
trajin test.mol2
box x 20 y 20 z 20 alpha 90 beta 90 gamma 90
energy out ewald.dat etype pme cut 5.6 dsumtol 0.0000001 skinnb 0.01
#vector UX ucellx
#vector UY ucelly
#vector UZ ucellz
#run
#writedata ucell.mol2 vectraj trajfmt mol2 UX UY UZ
EOF
    RunCpptraj "$UNITNAME"
  fi
}

NaCl() {
  UNITNAME='Particle mesh Ewald test (NaCl crystal)'
  CheckFor maxthreads 1
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm ../Test_Ewald/nacl.box.parm7
trajin ../Test_Ewald/nacl.box.rst7
debug actions 1
energy Reg nonbond out ewald.dat etype ewald cut 5.6 dsumtol 0.0000001 rsumtol 0.000000001 skinnb 0.01 mlimits 12,12,12

energy Pme nonbond out ewald.dat etype pme cut 5.6 dsumtol 0.0000001 skinnb 0.01 nfft 32,32,32
EOF
    RunCpptraj "$UNITNAME"
    grep "DEBUG: Eself" test.out > pme.nacl.dat
    DoTest pme.nacl.dat.save pme.nacl.dat
    DoTest ewald.dat.save ewald.dat
  fi
}

TrpzipNonortho() {
  UNITNAME='Particle mesh Ewald test (trunc. oct)'
  CheckFor netcdf maxthreads 1
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
debug actions 1
energy Reg out ew_tz2.dat etype ewald skinnb 0.01 \
       cut 8.0 dsumtol 0.0000001 rsumtol 0.000000001
energy Pme out ew_tz2.dat etype pme   skinnb 0.01 order 6 \
       cut 8.0 dsumtol 0.0000001 nfft 96,90,90
EOF
    RunCpptraj "$UNITNAME"
    #DoTest ew_tz2.dat.save ew_tz2.dat
  fi
}

TrpzipOrtho() {
  UNITNAME='Particle mesh Ewald test (ortho)'
  CheckFor netcdf maxthreads 1
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
debug actions 1
energy Reg nonbond out ew_tz2o.dat etype ewald skinnb 0.01 \
       cut 8.0 dsumtol 0.0000001 rsumtol 0.000000001
energy Pme nonbond out ew_tz2o.dat etype pme   skinnb 0.01 order 6 \
       cut 8.0 dsumtol 0.0000001 nfft 72,90,72
precision ew_tz2o.dat 20 10
EOF
    RunCpptraj "$UNITNAME"
    DoTest ew_tz2o.dat.save ew_tz2o.dat
  fi
}

Tz2_Nonortho_10() {
  UNITNAME='Ewald test (trunc. oct), 10 frames'
  CheckFor netcdf long
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
energy out ew_tz2_10.dat etype ewald skinnb 0.01
EOF
    RunCpptraj "$UNITNAME"
    DoTest ew_tz2_10.dat.save ew_tz2_10.dat
  fi
}

Tz2_Ortho_10() {
  UNITNAME='Ewald test (ortho), 10 frames'
  CheckFor netcdf long
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
energy out tz2_ortho.dat etype ewald skinnb 0.01
EOF
    RunCpptraj "Ewald test (ortho), 10 frames"
    DoTest tz2_ortho.dat.save tz2_ortho.dat
  fi
}

#Simple
NaCl
#TrpzipNonortho
TrpzipOrtho
#Tz2_Nonortho_10
#Tz2_Ortho_10

EndTest
exit 0 
