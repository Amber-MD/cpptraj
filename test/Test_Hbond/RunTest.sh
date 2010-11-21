#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles hbond.in nhb.dat avghb.dat 

# Test 1
CheckNetcdf
cat > hbond.in <<EOF
noprogress
parm ../DPDP.mod.GA12.parm7
trajin ../DPDP.nc
hbond out nhb.dat avgout avghb.dat 
EOF
INPUT="-i hbond.in"
RunCpptraj "Hbond command test."
DoTest nhb.dat.save nhb.dat 
DoTest avghb.dat.save avghb.dat

CheckTest

EndTest

exit 0
