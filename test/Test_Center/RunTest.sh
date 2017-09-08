#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles center.in centered.crd origin.centered.crd origin.mass.centered.crd \
           ref.centered.crd point.centered.crd
# NOTE: Also tests strip functionality
TESTNAME='Center tests'
Requires netcdf maxthreads 2

INPUT="-i center.in"
# Box center
cat > center.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
strip !(:1-23)
center :1-13
trajout centered.crd
EOF
RunCpptraj "Center command test."
DoTest centered.crd.save centered.crd
# Origin center
cat > center.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
strip !(:1-23)
center :1-13 origin
trajout origin.centered.crd
EOF
RunCpptraj "Center origin command test."
DoTest origin.centered.crd.save origin.centered.crd
# Origin center, mass-weighted
cat > center.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
strip !(:1-23)
center :1-13 origin mass
trajout origin.mass.centered.crd
EOF
RunCpptraj "Center origin mass command test."
DoTest origin.mass.centered.crd.save origin.mass.centered.crd
# Reference center
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
RunCpptraj "Center reference command test."
DoTest ref.centered.crd.save ref.centered.crd
# Center to point
cat > center.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 2
strip !:4
center point 2 4 6
trajout point.centered.crd
EOF
RunCpptraj "Center point test."
DoTest point.centered.crd.save point.centered.crd
EndTest
exit 0
