#!/bin/bash

. ../MasterTest.sh

CleanFiles box.in addbox.rst7 modX.rst7

INPUT="-i box.in"
cat > box.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
box x 42.428  y 42.428  z 42.428 alpha 109.471 beta 109.471 gamma 109.471
trajout addbox.rst7 restart novelocity time0 1.0
go
EOF
RunCpptraj "Box Test (Add box info)"
DoTest addbox.rst7.save addbox.rst7

cat > box.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
box x 45.0
trajout modX.rst7 restart time0 1.0
go
EOF
RunCpptraj "Box test (Modify box length)"
DoTest modX.rst7.save modX.rst7

#cat > box.in <<EOF
#parm ../tz2.parm7
#reference ../tz2.rst7
#trajin ../tz2.pdb 
#box x 42.428  y 42.428  z 42.428 alpha 90.0 beta 90.0 gamma 90.90
##image
#unwrap reference
#EOF
#RunCpptraj "Box test with image"

EndTest

exit 0
