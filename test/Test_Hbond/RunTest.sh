#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles hbond.in nhb.dat avghb.dat solvhb.dat solvavg.dat

# Solute test 
CheckNetcdf
cat > hbond.in <<EOF
noprogress
parm ../DPDP.parm7
trajin ../DPDP.nc
hbond out nhb.dat avgout avghb.dat 
EOF
INPUT="-i hbond.in"
RunCpptraj "Solute Hbond test."
DoTest nhb.dat.save nhb.dat 
DoTest avghb.dat.save avghb.dat

# Solvent test
cat > hbond.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
hbond hb out solvhb.dat avgout solvavg.dat :1-13 solventacceptor :WAT@O solventdonor :WAT \
      solvout solvavg.dat bridgeout solvavg.dat
EOF
RunCpptraj "Solvent Hbond test."
DoTest solvhb.dat.save solvhb.dat
DoTest solvavg.dat.save solvavg.dat

CheckTest

EndTest

exit 0
