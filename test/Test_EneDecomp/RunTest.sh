#!/bin/bash

. ../MasterTest.sh

TESTNAME='Energy decomposition tests'

CleanFiles enedecomp.in ene.*.dat decomp.*.dat Total.*.dat

TESTNAME='Particle mesh Ewald tests'
#Requires maxthreads 1 

INPUT='-i enedecomp.in'

UNITNAME='NaCl box decomposition'
CheckFor libpme maxthreads 1
if [ $? -eq 0 ] ; then
  cat > enedecomp.in <<EOF
parm ../Test_Ewald/nacl.box.parm7
trajin ../Test_Ewald/nacl.box.rst7
energy Pme nonbond out ene.nacl.box.dat \
  etype pme cut 5.6 dsumtol 0.0000001 skinnb 0.01 nfft 32,32,32

enedecomp ATM * out decomp.nacl.box.dat \
  pme cut 5.6 dsumtol 0.0000001 skinnb 0.01 nfft 32,32,32

#enedecomp ATM @4 out decomp.nacl.box.dat
run
Total = sum(ATM)
writedata Total.nacl.box.dat Total
EOF
  RunCpptraj "$UNITNAME"
  DoTest decomp.nacl.box.dat.save decomp.nacl.box.dat
  DoTest ene.nacl.box.dat.save ene.nacl.box.dat
  DoTest Total.nacl.box.dat.save Total.nacl.box.dat
fi

UNITNAME='AFV decomposition'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > enedecomp.in <<EOF
parm ../AFV.parm7 
trajin AFV.rst7
energy ENE out ene.AFV.dat
enedecomp ATM * out decomp.AFV.dat
#enedecomp ATM @4 out decomp.AFV.dat
run
Total = sum(ATM)
writedata Total.AFV.dat Total
EOF
  RunCpptraj "$UNITNAME"
  DoTest decomp.AFV.dat.save decomp.AFV.dat
  DoTest ene.AFV.dat.save ene.AFV.dat
  DoTest Total.AFV.dat.save Total.AFV.dat
fi

UNITNAME='Trpzip2 decomposition'
CheckFor netcdf maxthreads 10
if [ $? -eq 0 ] ; then
  cat > enedecomp.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
enedecomp TZ2 * out decomp.tz2.dat
run
EOF
  RunCpptraj "$UNITNAME"
  DoTest decomp.tz2.dat.save decomp.tz2.dat
fi

EndTest
