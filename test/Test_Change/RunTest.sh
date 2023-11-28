#!/bin/bash

. ../MasterTest.sh

CleanFiles change.in ala3.mod.pdb ala3.chain.pdb crdala3.chain.pdb \
           AFV.zeroHmass.dat AFV.fluctMass.dat merged.?-?.mol2

TESTNAME='Change command test'

INPUT='-i change.in'

UNITNAME='Change atom and residue name tests'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > change.in <<EOF
parm ../Test_Charmm/ala3.psf
trajin ../Test_Charmm/ala3.dcd 1 1
change parmindex 0 resname from :1 to NALA
change parmindex 0 resname from :3 to CALA
change parmindex 0 atomname from @HN to H
trajout ala3.mod.pdb chainid " "
EOF
  RunCpptraj "$UNITNAME"
  DoTest ala3.mod.pdb.save ala3.mod.pdb
fi

UNITNAME='Change chain ID test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > change.in <<EOF
parm ../Test_Charmm/ala3.psf
trajin ../Test_Charmm/ala3.dcd 1 1
change parmindex 0 chainid of * to A
trajout ala3.chain.pdb
EOF
  RunCpptraj "$UNITNAME"
  DoTest ala3.chain.pdb.save ala3.chain.pdb
fi

UNITNAME='Change chain ID of COORDS set test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > change.in <<EOF
parm ../Test_Charmm/ala3.psf
loadcrd ../Test_Charmm/ala3.dcd 1 1 name MyCrd
change crdset MyCrd chainid of * to A
crdout MyCrd crdala3.chain.pdb
EOF
  RunCpptraj "$UNITNAME"
  DoTest ala3.chain.pdb.save crdala3.chain.pdb
fi

UNITNAME='Change mass of atoms test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > change.in <<EOF
parm ../AFV.parm7
change mass of @H= to 0.0
atoms * out AFV.zeroHmass.dat
#parmwrite out AFV.zeroHmass.parm7
EOF
  RunCpptraj "$UNITNAME"
  DoTest AFV.zeroHmass.dat.save AFV.zeroHmass.dat
fi

UNITNAME='Change mass from data set test'
CheckFor maxthreads 1 # FIXME this should eventually work in parallel
if [ $? -eq 0 ] ; then
  cat > change.in <<EOF
parm ../AFV.parm7
trajin ../AFV.nc
align first !@H=
atomicfluct MyFluct #out AFV.myfluct.dat
run
change mass fromset MyFluct
atoms * out AFV.fluctMass.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest AFV.fluctMass.dat.save AFV.fluctMass.dat
fi

UNITNAME='Merge residues test'
cat > change.in <<EOF
set DIR = ../Test_Sequence
parm \$DIR/Nucleotide.ic.charge.mol2.save

loadcrd \$DIR/Nucleotide.ic.charge.mol2.save name MyCrd1
loadcrd \$DIR/Nucleotide.ic.charge.mol2.save name MyCrd2
loadcrd \$DIR/Nucleotide.ic.charge.mol2.save name MyCrd3

change crdset MyCrd1 mergeres firstres 1 lastres 3
change crdset MyCrd2 mergeres firstres 1 lastres 2
change crdset MyCrd3 mergeres firstres 2 lastres 3

crdout MyCrd1 merged.1-3.mol2
crdout MyCrd2 merged.1-2.mol2
crdout MyCrd3 merged.2-3.mol2
EOF
RunCpptraj "$UNITNAME"
DoTest merged.1-3.mol2.save merged.1-3.mol2
DoTest merged.1-2.mol2.save merged.1-2.mol2
DoTest merged.2-3.mol2.save merged.2-3.mol2

EndTest
exit 0
