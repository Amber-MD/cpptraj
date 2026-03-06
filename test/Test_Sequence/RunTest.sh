#!/bin/bash

. ../MasterTest.sh

TESTNAME='Sequence tests'

CleanFiles cpptraj.in Mol.mol2 Mol2.mol2 MOC.mol2 CNALA.mol2 \
           Mol3.mol2 Mol4.mol2 Nucleotide.ic.charge.mol2

INPUT='-i cpptraj.in'

Basic() {
  cat > cpptraj.in <<EOF
readdata ../Test_ReadOFF/aminocn15ipq_10.0.lib name A15
readdata cph_nucleic_caps.lib name CAPS
list
#sequence CAPS[MOC] A15[CNALA] name Mol
#crdout Mol Mol.mol2

sequence MOC CNALA name Mol2
crdout Mol2 Mol2.mol2

crdout CAPS[MOC] MOC.mol2
crdout A15[CNALA] CNALA.mol2
EOF
  RunCpptraj "$TESTNAME, library files"
#  DoTest Mol.mol2.save Mol.mol2
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
}

# Link nucleic acid base + sugar + phosphate, IC, fix charges.
# NOTE: This is not really how something like this should be
#       constructed since after atoms are stripped from the sugar
#       the original chirality is lost; sequence then builds in a
#       default chirality which is not correct for nucleic acids.
#       This test is purely a regression test.
DNAic_charge() {
  UNITNAME="$TESTNAME, Construct fake Nucleic Acid, IC, fix charges"
  CheckFor maxthreads 1
  if [ $? -eq 0 ] ; then
    cat > cpptraj.in <<EOF
set DIR = ../Test_Graft
parm \$DIR/DDD.names.mol2
loadcrd \$DIR/DDD.names.mol2 name Sugar parm DDD.names.mol2
parm \$DIR/MP1.names.mol2
loadcrd \$DIR/MP1.names.mol2 name Phos parm MP1.names.mol2
parm \$DIR/ADD.names.mol2
loadcrd \$DIR/ADD.names.mol2 name Base parm ADD.names.mol2
#list
# Strip components
crdaction Base  strip charge  0.116  @C1,H1,H6,H7
crdaction Sugar keep  charge  0.2933 keepmask !(@C4,H6,H7,H8,H12,H11,O1,H1)
crdaction Phos  strip charge -1.4093 @C3,H4,H5,H6,O1,C1,H1,H2,H3
# Set connect atoms
dataset connect Base               tailmask @N1
dataset connect Sugar headmask @C3 tailmask @C1
dataset connect Phos  headmask @O3
sequence name MyMol Base Sugar Phos
charge crdset MyMol *
crdout MyMol Nucleotide.ic.charge.mol2
EOF
    RunCpptraj "$UNITNAME"
    #DoTest Nucleotide.ic.charge.mol2.save Nucleotide.ic.charge.mol2
    DoTest Nucleotide.wrongChirality.ic.charge.mol2.save Nucleotide.ic.charge.mol2
  fi
}

Basic
DNAic_charge

EndTest
