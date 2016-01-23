#!/bin/bash

. ../MasterTest.sh

CleanFiles les.in splittraj.? avg.crd
MaxThreads 5 "LES split/avg tests"
if [[ $? -ne 0 ]] ; then
  EndTest
  exit 0
fi
INPUT="-i les.in"
cat > les.in <<EOF
parm lesparm
trajin run0.crd
lessplit out splittraj average avg.crd
EOF
RunCpptraj "LES split test"
DoTest splittraj.0.save splittraj.0
DoTest avg.crd.save avg.crd
EndTest
exit 0
