#!/bin/bash

. ../MasterTest.sh

CleanFiles les.in splittraj.? avg.crd
TESTNAME='LES split/average test'
Requires maxthreads 5

INPUT="-i les.in"
cat > les.in <<EOF
parm lesparm
trajin run0.crd
lessplit out splittraj average avg.crd
EOF
RunCpptraj "$TESTNAME"
DoTest splittraj.0.save splittraj.0
DoTest avg.crd.save avg.crd
EndTest
exit 0
