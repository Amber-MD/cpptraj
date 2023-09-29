#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in zmatrix.dat

INPUT='-i cpptraj.in'

TESTNAME='Zmatrix tests'

cat > cpptraj.in <<EOF
parm ../Test_ReadPrep/test1.mol2.save
loadcrd ../Test_ReadPrep/test1.mol2.save name MyCrd
zmatrix MyCrd out zmatrix.dat
EOF
RunCpptraj "$TESTNAME"
DoTest zmatrix.dat.save zmatrix.dat

EndTest
