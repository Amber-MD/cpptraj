#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in tz2.hmr.parm7

TESTNAME='Hydrogen mass repartition test'

INPUT='-i cpptraj.in'
cat > cpptraj.in <<EOF
parm ../tz2.parm7
hmassrepartition
parmwrite out tz2.hmr.parm7
EOF
RunCpptraj "$TESTNAME"
DoTest tz2.hmr.parm7.save tz2.hmr.parm7 -I %VERSION

EndTest
