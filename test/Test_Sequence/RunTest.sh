#!/bin/bash

. ../MasterTest.sh

TESTNAME='Sequence tests'

CleanFiles cpptraj.in Mol.mol2 Mol2.mol2 MOC.mol2 CNALA.mol2 \
           Mol3.mol2 Mol4.mol2

INPUT='-i cpptraj.in'

cat > cpptraj.in <<EOF
readdata ../Test_ReadOFF/aminocn15ipq_10.0.lib name A15
readdata cph_nucleic_caps.lib name CAPS
list
sequence CAPS[MOC] A15[CNALA] name Mol
crdout Mol Mol.mol2

sequence libset A15 libset CAPS MOC CNALA name Mol2
crdout Mol2 Mol2.mol2

crdout CAPS[MOC] MOC.mol2
crdout A15[CNALA] CNALA.mol2
EOF
RunCpptraj "$TESTNAME, library files"
DoTest Mol.mol2.save Mol.mol2
DoTest Mol.mol2.save Mol2.mol2

# NOTE: Depends on mol2 generation of previous test
cat > cpptraj.in <<EOF
parm MOC.mol2
loadcrd MOC.mol2 parm MOC.mol2 name MOC
dataset connect MOC tail 5
parm CNALA.mol2
loadcrd CNALA.mol2 parm CNALA.mol2 name CNALA
dataset connect CNALA head 1
sequence MOC CNALA name Mol3
crdout Mol3 Mol3.mol2

dataset connect MOC tailmask @O5
dataset connect CNALA headmask @N
sequence MOC CNALA name Mol4
crdout Mol4 Mol4.mol2
EOF
RunCpptraj "$TESTNAME, mol2 files"
DoTest Mol.mol2.save Mol3.mol2
DoTest Mol.mol2.save Mol4.mol2

EndTest
