#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles mask.in mask.out mask.pdb.1 mask.mol2.1 M.dat

CheckNetcdf

INPUT="-i mask.in"
# Test 1
MaxThreads 1 "Mask single frame tests"
if [[ $? -eq 0 ]] ; then
  cat > mask.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
mask "(:5 <:3.0) & :WAT" maskout mask.out maskpdb mask.pdb
mask "(:5 <:3.0) & :WAT" maskmol2 mask.mol2
EOF
RunCpptraj "Mask command test."
DoTest mask.out.save mask.out
DoTest mask.pdb.1.save mask.pdb.1
DoTest mask.mol2.1.save mask.mol2.1
fi

cat > mask.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
mask "(:8@NZ <:3.0) & :WAT@O" name M out M.dat noxcol
EOF
RunCpptraj "Mask longer command test."
DoTest M.dat.save M.dat

EndTest

exit 0
