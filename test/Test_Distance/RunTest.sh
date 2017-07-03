#!/bin/bash

. ../MasterTest.sh

CleanFiles dist.in dist.dat

INPUT='-i dist.in'

cat > dist.in <<EOF
parm ../tz2.parm7
reference ../tz2.pdb
trajin ../tz2.nc

distance ToRef @1 @1 out dist.dat reference
EOF
RunCpptraj "Distance tests."
DoTest dist.dat.save dist.dat

EndTest
exit 0
