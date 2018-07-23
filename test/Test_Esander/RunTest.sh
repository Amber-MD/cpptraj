#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in Esander.dat force.nc Edpdp.dat CpptrajEsander.parm7 NoWat.dat
TESTNAME='SANDER energy tests (long)'
Requires sanderlib netcdf long

INPUT="-i ene.in"

TestPME() {
  UNITNAME='SANDER energy test, PME'
  CheckFor pnetcdf maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
esander S out Esander.dat saveforces
trajout force.nc
EOF
    RunCpptraj "$UNITNAME"
    DoTest Esander.dat.save Esander.dat
    NcTest force.nc.save force.nc -a 0.000001
  fi
}

TestGB() {
  cat > ene.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
esander DPDP out Edpdp.dat gbsa 1
EOF
  RunCpptraj "SANDER energy test, GB."
  DoTest Edpdp.dat.save Edpdp.dat
}

TestStrip() {
  UNITNAME="SANDER energy test after 'strip', PME."
  CheckFor maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
strip :WAT
esander NoWat out NoWat.dat
EOF
    RunCpptraj "$UNITNAME"
    DoTest NoWat.dat.save NoWat.dat
  fi
}

TestGB
TestPME
TestStrip

EndTest
exit 0
