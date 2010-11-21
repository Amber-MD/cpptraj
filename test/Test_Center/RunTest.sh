#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles center.in centered.crd origin.centered.crd origin.mass.centered.crd

# Test 1
CheckNetcdf
cat > center.in <<EOF
noprogress
parm ../ChainA-tip3p.parm7
trajin ../run0.nc 1 10
center :1-268
trajout centered.crd
EOF
INPUT="-i center.in"
RunCpptraj "Center command test."
DoTest centered.crd.save centered.crd

cat > center.in <<EOF
noprogress
parm ../ChainA-tip3p.parm7
trajin ../run0.nc 1 10
center :1-268 origin
trajout origin.centered.crd
EOF
INPUT="-i center.in"
RunCpptraj "Center origin command test."
DoTest origin.centered.crd.save origin.centered.crd

cat > center.in <<EOF
noprogress
parm ../ChainA-tip3p.parm7
trajin ../run0.nc 1 10
center :1-268 origin mass
trajout origin.mass.centered.crd
EOF
INPUT="-i center.in"
RunCpptraj "Center origin mass command test."
DoTest origin.mass.centered.crd.save origin.mass.centered.crd

CheckTest

EndTest

exit 0
