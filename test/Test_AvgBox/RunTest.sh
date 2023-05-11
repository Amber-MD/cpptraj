#!/bin/bash

. ../MasterTest.sh

CleanFiles box.in box.dat

INPUT='-i box.in'

TESTNAME='Average box tests'
Requires netcdf

cat > box.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
avgbox MyBox out box.dat
EOF
RunCpptraj "Average box test"
DoTest box.dat.save box.dat

EndTest
exit 0
