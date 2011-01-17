#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles closest.in Closest.pdb closest2.in first.Closest.pdb 

# Test 1
CheckNetcdf
cat > closest.in <<EOF
noprogress
parm ../ChainA-tip3p.parm7
trajin ../run0.nc 1 1
closest 10 :1-268 first
trajout first.Closest.pdb pdb
EOF
INPUT="-i closest.in"
RunCpptraj "Closest command test using first solvent atom."
DoTest first.Closest.pdb.save first.Closest.pdb

# Long test, disabled for now.
#cat > closest2.in <<EOF
#noprogress
#parm ../ChainA-tip3p.parm7
#trajin ../run0.nc 1 1
#closest 10 :1-268
#trajout Closest.pdb pdb
#EOF
#INPUT="-i closest2.in"
#RunCpptraj "Closest command test using all solvent atoms."
#DoTest Closest.pdb.save Closest.pdb

CheckTest

EndTest

exit 0
