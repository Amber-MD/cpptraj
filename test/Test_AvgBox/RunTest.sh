#!/bin/bash

. ../MasterTest.sh

CleanFiles box.in box.ortho.dat

INPUT='-i box.in'

TESTNAME='Average box tests'
Requires netcdf

cat > box.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
avgbox MyBox out box.ortho.dat
EOF
RunCpptraj "Average box test"
DoTest box.ortho.dat.save box.ortho.dat

EndTest
exit 0
