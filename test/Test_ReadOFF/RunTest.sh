#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in Lib.mol2

INPUT='-i cpptraj.in'

cat > cpptraj.in <<EOF
readdata aminocn15ipq_10.0.lib name Lib
crdout Lib[CNALA] Lib.mol2
list
quit
EOF
RunCpptraj "Amber OFF format read test"
DoTest Lib.mol2.save Lib.mol2

EndTest
