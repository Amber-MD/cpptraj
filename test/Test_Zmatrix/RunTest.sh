#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in zmatrix.dat zmatrix.0SB.dat

INPUT='-i cpptraj.in'

TESTNAME='Zmatrix tests'

Basic() {
  cat > cpptraj.in <<EOF
parm ../Test_ReadPrep/test1.mol2.save
loadcrd ../Test_ReadPrep/test1.mol2.save name MyCrd
zmatrix MyCrd out zmatrix.dat
EOF
  RunCpptraj "$TESTNAME, Basic"
  DoTest zmatrix.dat.save zmatrix.dat
}

Sugar() {
  cat > cpptraj.in <<EOF
parm 0SB.mol2
loadcrd 0SB.mol2 name 0SB
zmatrix 0SB out zmatrix.0SB.dat
EOF
  RunCpptraj "$TESTNAME, Sugar"
  DoTest zmatrix.0SB.dat.save zmatrix.0SB.dat
}

Basic
Sugar

EndTest
