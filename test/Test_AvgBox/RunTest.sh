#!/bin/bash

. ../MasterTest.sh

CleanFiles box.in box.ortho.dat box.truncoct.dat

INPUT='-i box.in'

TESTNAME='Average box tests'
Requires netcdf maxthreads 10

cat > box.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
avgbox MyBox out box.ortho.dat
EOF
RunCpptraj "Average box test, orthogonal"
DoTest box.ortho.dat.save box.ortho.dat

cat > box.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
avgbox MyBox out box.truncoct.dat
EOF
RunCpptraj "Average box test, truncated octahedron"
DoTest box.truncoct.dat.save box.truncoct.dat

EndTest
exit 0
