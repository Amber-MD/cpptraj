#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in Esander.dat force.nc

INPUT="-i ene.in"

TestPME() {
  cat > ene.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
esander S out Esander.dat saveforces
trajout force.nc
EOF
  RunCpptraj "Libsander test."
  DoTest Esander.dat.save Esander.dat
  NcTest force.nc.save force.nc
}

TestPME

EndTest
exit 0
