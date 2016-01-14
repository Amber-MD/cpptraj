#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles center.in centered.crd origin.centered.crd origin.mass.centered.crd ref.centered.crd
# Also tests strip functionality
MaxThreads 2 "Center tests."
if [[ $? -ne 0 ]] ; then
  echo ""
  exit 0
fi

# Test 1
CheckNetcdf
cat > center.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
strip !(:1-23)
center :1-13
trajout centered.crd
EOF
INPUT="-i center.in"
RunCpptraj "Center command test."
DoTest centered.crd.save centered.crd

cat > center.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
strip !(:1-23)
center :1-13 origin
trajout origin.centered.crd
EOF
INPUT="-i center.in"
RunCpptraj "Center origin command test."
DoTest origin.centered.crd.save origin.centered.crd

cat > center.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
strip !(:1-23)
center :1-13 origin mass
trajout origin.mass.centered.crd
EOF
INPUT="-i center.in"
RunCpptraj "Center origin mass command test."
DoTest origin.mass.centered.crd.save origin.mass.centered.crd

cat > center.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
strip !(:1-23)
parm ../tz2.parm7 [NOWAT]
reference ../tz2.nc 1 parm [NOWAT]
center :1-13 reference :1-12
trajout ref.centered.crd
EOF
INPUT="-i center.in"
RunCpptraj "Center reference command test."
DoTest ref.centered.crd.save ref.centered.crd
CheckTest

EndTest

exit 0
