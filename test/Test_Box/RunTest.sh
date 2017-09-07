#!/bin/bash

. ../MasterTest.sh

CleanFiles box.in addbox.rst7 addbox.rst7.? addbox.rst7.10 \
                  modX.rst7   modX.rst7.?   modX.rst7.10
RequiresNetcdf "Box tests"
RequiresMaxThreads 10 "Box tests"

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

EndTest

exit 0
