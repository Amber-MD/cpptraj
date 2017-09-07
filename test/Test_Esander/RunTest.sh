#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in Esander.dat force.nc Edpdp.dat CpptrajEsander.parm7 NoWat.dat

RequiresSanderlib "SANDER energy tests"
RequiresNetcdf "SANDER energy tests"

INPUT="-i ene.in"

TestPME() {
  CheckPnetcdf "SANDER energy test, PME"
  if [[ $? -eq 0 ]] ; then
    cat > ene.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
esander S out Esander.dat saveforces
trajout force.nc
EOF
    RunCpptraj "SANDER energy test, PME."
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
  cat > ene.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
strip :WAT
esander NoWat out NoWat.dat
EOF
  RunCpptraj "SANDER energy test after 'strip', PME."
  DoTest NoWat.dat.save NoWat.dat
}

TestGB
TestPME
TestStrip

EndTest
exit 0
