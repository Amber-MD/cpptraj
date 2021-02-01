#!/bin/bash

. ../MasterTest.sh

CleanFiles dist.in dist.dat ortho.dat truncoct.dat

TESTNAME='Distance tests'
Requires netcdf
INPUT='-i dist.in'

cat > dist.in <<EOF
parm ../tz2.parm7
reference ../tz2.pdb
trajin ../tz2.nc

distance EndToEnd :1 :13 out dist.dat
distance ToRef @1 @1 out dist.dat reference
distance Point :1 point 0.0 0.0 0.0 out dist.dat
run
EOF
RunCpptraj "$TESTNAME"
DoTest dist.dat.save dist.dat

UNITNAME='Orthogonal imaged distance'
cat > dist.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
distance EndToEnd :1 :13 out ortho.dat
EOF
RunCpptraj "$UNITNAME"
DoTest ortho.dat.save ortho.dat

UNITNAME='Non-orthogonal imaged distance'
cat > dist.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
distance EndToEnd :1 :13 out truncoct.dat
EOF
RunCpptraj "$UNITNAME"
DoTest truncoct.dat.save truncoct.dat

EndTest
exit 0
