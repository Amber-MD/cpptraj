#!/bin/bash

. ../MasterTest.sh

TESTNAME='XYZ format tests'

CleanFiles xyz.in tz2.xyz test1.crd.save test1.crd

INPUT='-i xyz.in'

UNITNAME='XYZ format write'
cat > xyz.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd 1 10
trajout test1.crd.save crd
trajout tz2.xyz
EOF
RunCpptraj "$UNITNAME"

UNITNAME='Atom XYZ format read'
cat > xyz.in <<EOF
parm ../tz2.parm7
trajin tz2.xyz
trajout test1.crd
EOF
RunCpptraj "$UNITNAME"
DoTest test1.crd.save test1.crd

EndTest
exit 0
