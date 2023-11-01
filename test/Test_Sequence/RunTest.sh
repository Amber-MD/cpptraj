#!/bin/bash

. ../MasterTest.sh

TESTNAME='Sequence tests'

CleanFiles cpptraj.in Mol.mol2

INPUT='-i cpptraj.in'

cat > cpptraj.in <<EOF
readdata ../Test_ReadOFF/aminocn15ipq_10.0.lib name A15
readdata cph_nucleic_caps.lib name CAPS
list
sequence CAPS[MOC] A15[CNALA] name Mol
crdout Mol Mol.mol2
EOF
RunCpptraj "$TESTNAME"
DoTest Mol.mol2.save Mol.mol2

EndTest
