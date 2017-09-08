#!/bin/bash

. ../MasterTest.sh

CleanFiles dist.in dist.dat

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

EndTest
exit 0
