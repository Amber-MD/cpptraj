#!/bin/bash

. ../MasterTest.sh

CleanFiles esander.in gb.dat pme.dat Esander.dat force.crd Edpdp.dat CpptrajEsander.parm7 NoWat.dat

INPUT='-i esander.in'
TESTNAME='Energy via sander API short tests (GB/PME)'
Requires sanderlib maxthreads 10

UNITNAME='Short libsander tests'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > esander.in <<EOF
parm ../Test_Hbond/strip.4lztSc_nowat.parm7
trajin ../Test_Hbond/strip.4lztSc.rst7
esander 4lztSc out gb.dat \
  saltcon 1.0 cut 9999.0 gbsa 1 igb 1 ntb 0
esander PME out pme.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest gb.dat.save gb.dat
  DoTest pme.dat.save pme.dat
fi

# Set to 1 if compiled with -DCPPTRAJ_ESANDER_ENERGY_ORDER
NEWORDER=0
if [ $NEWORDER -eq 1 ] ; then
  SAVEPREFIX='neworder.'
else
  SAVEPREFIX=''
fi

TestPME() {
  UNITNAME='Short SANDER energy test, PME'
  CheckFor netcdf maxthreads 3
  if [ $? -eq 0 ] ; then
    cat > esander.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 3
esander S out Esander.dat saveforces
trajout force.crd
EOF
    RunCpptraj "$UNITNAME"
    DoTest "$SAVEPREFIX"Esander.dat.save Esander.dat
    DoTest force.crd.save force.crd
  fi
}

TestGB() {
  UNITNAME="Short SANDER energy test, GB"
  CheckFor netcdf maxthreads 10
  cat > esander.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc 1 10
esander DPDP out Edpdp.dat gbsa 1
EOF
  RunCpptraj "$UNITNAME"
  DoTest "$SAVEPREFIX"Edpdp.dat.save Edpdp.dat
}

TestGB
TestPME

EndTest
exit 0
