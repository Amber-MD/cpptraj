#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in temperature.dat rmsd.dat ene.dat offset.*.dat

TESTNAME='TNG read test'

Requires tng maxthreads 10

INPUT='-i cpptraj.in'

TestTng() {
  UNITNAME=$TESTNAME
  trajin_args=$1
  refcmd=''
  rmsref='first'
  xmin=0
  xstep=10
  tout='temperature.dat'
  rout='rmsd.dat'
  eout='ene.dat'
  skiptest=0
  if [ ! -z "$trajin_args" ] ; then
    UNITNAME="$TESTNAME (with offset)"
    refcmd='reference md_1_1.tng 1'
    rmsref='reference'
    xmin=10
    xstep=20
    tout='offset.temperature.dat'
    rout='offset.rmsd.dat'
    eout='offset.ene.dat'
    CheckFor maxthreads 4
    skiptest=$?
  fi
  if [ $skiptest -eq 0 ] ; then
    cat > cpptraj.in <<EOF
parm topol.parm7
noprogress
$refcmd
trajin md_1_1.tng $trajin_args
# Temperature calc tests velocity read
temperature MyTemp ntc 2 out $tout xmin $xmin xstep $xstep prec 12.7
# Energy (specifically kinetic VV) tests velocity/force read.
# It is probably not correct since this really requires the half
# step velocities, but is good enough for a regression test.
energy MyEne ^1 out $eout xmin $xmin xstep $xstep prec 12.7 \
  kinetic ketype vv dt 0.002 
# RMSD tests coordinates read
rms MyRms $rmsref ^1&@C,CA,N
run
# Convert RMSD to nm
MyRmsNm = MyRms / 10.0
writedata $rout xmin $xmin xstep $xstep MyRmsNm prec 12.7
EOF
    RunCpptraj "$UNITNAME"
    if [ -z "$trajin_args" ] ; then
      DoTest temperature.dat.save temperature.dat
      DoTest rmsd.dat.save rmsd.dat
      DoTest ene.dat.save ene.dat
    else
      DoTest offset.temperature.dat.save offset.temperature.dat
      DoTest offset.rmsd.dat.save offset.rmsd.dat
      DoTest offset.ene.dat.save offset.ene.dat
   fi
  fi
}

TestTng

TestTng "2 8 2"

EndTest

exit 0
