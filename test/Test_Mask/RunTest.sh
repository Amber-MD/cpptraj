#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles mask.in mask.out mask.pdb.1 mask.mol2.1

# Test 1
CheckNetcdf
cat > mask.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
mask "(:5 <:3.0) & :WAT" maskout mask.out maskpdb mask.pdb
mask "(:5 <:3.0) & :WAT" maskmol2 mask.mol2
EOF
INPUT="-i mask.in"
RunCpptraj "Mask command test."
DoTest mask.out.save mask.out
DoTest mask.pdb.1.save mask.pdb.1
DoTest mask.mol2.1.save mask.mol2.1

CheckTest

EndTest

exit 0
