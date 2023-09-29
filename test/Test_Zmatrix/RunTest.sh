#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in zmatrix.dat zmatrix.0SB.dat fromZmatrix.0SB.mol2 \
           zmatrix.MEX.dat fromZmatrix.MEX.mol2

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
zmatrix 0SB out zmatrix.0SB.dat name MyZ
zmatrix 0SB zset MyZ name MyCrd
crdout MyCrd fromZmatrix.0SB.mol2
EOF
  RunCpptraj "$TESTNAME, Sugar"
  DoTest zmatrix.0SB.dat.save zmatrix.0SB.dat
  DoTest 0SB.mol2 fromZmatrix.0SB.mol2
}

MEX() {
  cat > cpptraj.in <<EOF
parm MEX.mol2
loadcrd MEX.mol2 name MEX
zmatrix MEX out zmatrix.MEX.dat name MyZ
zmatrix MEX zset MyZ name MyCrd
crdout MyCrd fromZmatrix.MEX.mol2
EOF
  RunCpptraj "$TESTNAME, methyl"
  DoTest MEX.mol2 fromZmatrix.MEX.mol2
}

Basic
Sugar
MEX

EndTest
