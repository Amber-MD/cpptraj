#!/bin/bash

. ../MasterTest.sh

CleanFiles translate.in translate.2.11.mol2 trp.topoint.1.2.3.mol2

TESTNAME='Translate test'
Requires maxthreads 1

INPUT='-i translate.in'

UNITNAME='Translate coordinates by delta test'
cat > translate.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
strip !(:2,11)
translate :2 x  1.0 y  10.0 z -0.5
translate :1 x  2.0 y -10.0 z -0.5
trajout translate.2.11.mol2
EOF
RunCpptraj "$UNITNAME"
DoTest translate.2.11.mol2.save translate.2.11.mol2

UNITNAME='Translate to specific coordinates test'
cat > translate.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
strip !:2
translate :1 topoint 1,2,3
trajout trp.topoint.1.2.3.mol2
EOF
RunCpptraj "$UNITNAME"
DoTest trp.topoint.1.2.3.mol2.save trp.topoint.1.2.3.mol2

EndTest
exit 0
