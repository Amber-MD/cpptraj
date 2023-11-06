#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in test1.mol2

INPUT='-i cpptraj.in'

TESTNAME='Test Prep File Read and Zmatrix Functionality'

cat > cpptraj.in <<EOF
readdata epACE.prepin name ACE
crdout ACE[*] test1.mol2
EOF
RunCpptraj "$TESTNAME"
DoTest test1.mol2.save test1.mol2

EndTest
