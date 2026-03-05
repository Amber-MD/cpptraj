#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in Esander.dat force.crd Edpdp.dat CpptrajEsander.parm7 NoWat.dat
TESTNAME='SANDER energy tests (short)'
Requires sanderlib netcdf

# Set to 1 if compiled with -DCPPTRAJ_ESANDER_ENERGY_ORDER
NEWORDER=0
if [ $NEWORDER -eq 1 ] ; then
  SAVEPREFIX='neworder.'
else
  SAVEPREFIX=''
fi

INPUT="-i ene.in"

TestPME() {
  UNITNAME='Short SANDER energy test, PME'
  CheckFor maxthreads 3
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 3
esander S out Esander.dat saveforces
trajout force.crd
EOF
    RunCpptraj "$UNITNAME"
    #DoTest "$SAVEPREFIX"Esander.dat.save Esander.dat
    #NcTest force.nc.save force.nc -a 0.000001
  fi
}

TestGB() {
  UNINTNAME="Short SANDER energy test, GB"
  CheckFor maxthreads 10
  cat > ene.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc 1 10
esander DPDP out Edpdp.dat gbsa 1
EOF
  RunCpptraj "$UNITNAME"
  #DoTest "$SAVEPREFIX"Edpdp.dat.save Edpdp.dat
}


TestGB
TestPME

EndTest
exit 0
