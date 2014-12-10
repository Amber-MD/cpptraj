#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles average.in test.pdb

CheckNetcdf
TOP="../tz2.parm7"

# Test 1
cat > average.in <<EOF
noprogress
trajin ../tz2.nc
average test.pdb pdb chainid X
EOF
INPUT="average.in"
RunCpptraj "Average Test."
DoTest test.pdb.save test.pdb

CheckTest

EndTest

exit 0
