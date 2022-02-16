#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles hbond.in nhb.dat avghb.dat solvhb.dat solvavg.dat bridges.dat

INPUT="-i hbond.in"

# Solute-solute, split by parts
TestUUsplit() {
  UNITNAME='Solute Hbond test with split'
  CheckFor netcdf
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
noprogress
parm ../DPDP.parm7
trajin ../DPDP.nc
hbond HB out nhb.dat avgout avghb.dat splitframe 50
run
EOF
    RunCpptraj "$UNITNAME"
    DoTest avghb.dat.save avghb.dat
  fi
}

# Solute-Solvent test, split by parts
TestUVsplit() {
  UNITNAME='Solute-solvent hbond test with split'
  CheckFor netcdf maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
hbond hb out solvhb.dat :1-13 solventacceptor :WAT@O solventdonor :WAT \
      solvout solvavg.dat bridgeout solvavg.dat splitframe 5 \
      bseries bseriesfile bridges.dat
run
EOF
    RunCpptraj "$UNITNAME"
    DoTest solvavg.dat.save solvavg.dat
    DoTest bridges.dat.save bridges.dat
  fi
}

TestUUsplit
TestUVsplit

EndTest

exit $? 
