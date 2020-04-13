#!/bin/bash

. ../MasterTest.sh

CleanFiles box.in addbox.rst7 addbox.rst7.? addbox.rst7.10 \
                  modX.rst7   modX.rst7.?   modX.rst7.10 \
                  frame1.rst7 tz2.box.rst7

TESTNAME='Box tests'
Requires netcdf maxthreads 10

INPUT="-i box.in"
cat > box.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10 
strip !(:1)
box x 42.428  y 42.428  z 42.428 alpha 109.471 beta 109.471 gamma 109.471
trajout addbox.rst7
go
EOF
RunCpptraj "Box Test (Add box info)"
DoTest addbox.rst7.1.save addbox.rst7.1
DoTest addbox.rst7.10.save addbox.rst7.10

cat > box.in <<EOF
parm ../tz2.parm7
parmstrip !(:1)
trajin addbox.rst7.1.save
box nobox
trajout frame1.rst7
EOF
RunCpptraj "Box test (remove box info)"
DoTest frame1.rst7.save frame1.rst7

cat > box.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 10
strip !(:1)
box x 45.0
trajout modX.rst7
go
EOF
RunCpptraj "Box test (Modify box length)"
DoTest modX.rst7.1.save modX.rst7.1
DoTest modX.rst7.10.save modX.rst7.10

cat > box.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 1
box auto 3.0
trajout tz2.box.rst7
run
EOF
RunCpptraj "Box test (auto orthogonal box)"
DoTest tz2.box.rst7.save tz2.box.rst7

EndTest

exit 0
