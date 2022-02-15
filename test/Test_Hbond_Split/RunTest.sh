#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles hbond.in nhb.dat avghb.dat

INPUT="-i hbond.in"

# Solute-solute, split by parts
TestUUsplit() {
  UNITNAME='Solute Hbond test'
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


TestUUsplit

EndTest

exit $? 
