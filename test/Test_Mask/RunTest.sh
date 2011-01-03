#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles mask.in mask.out 

# Test 1
CheckNetcdf
cat > mask.in <<EOF
noprogress
parm ../ChainA-tip3p.parm7
trajin ../run0.nc 1 1
mask "(:195 <:3.0) & :WAT" out mask.out
EOF
INPUT="-i mask.in"
RunCpptraj "Mask command test."
DoTest mask.out.save mask.out

CheckTest

EndTest

exit 0
