#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in Esander.dat force.nc Edpdp.dat

CheckSanderlib "SANDER energy tests"
if [[ $? -ne 0 ]] ; then
  EndTest
  exit 0
fi

INPUT="-i ene.in"

TestPME() {
  cat > ene.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
esander S out Esander.dat saveforces
trajout force.nc
EOF
  RunCpptraj "SANDER energy test, PME."
  DoTest Esander.dat.save Esander.dat
  NcTest force.nc.save force.nc
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

TestGB
TestPME

EndTest
exit 0
