#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles average.in test.pdb test2.pdb

RequiresNetcdf "Average test."
TOP="../tz2.parm7"

# Test 1
cat > average.in <<EOF
trajin ../tz2.nc
average test.pdb pdb chainid X
average crdset Tz2Avg 
run
crdout Tz2Avg test2.pdb chainid X
EOF
INPUT="average.in"
RunCpptraj "Average Test."
DoTest test.pdb.save test.pdb
DoTest test.pdb.save test2.pdb

CheckTest

EndTest

exit 0
