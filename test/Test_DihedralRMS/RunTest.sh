#!/bin/bash

. ../MasterTest.sh

CleanFiles dih.in 

TESTNAME='Dihedral RMSD tests'
Requires maxthreads 10

INPUT="-i dih.in"
UNITNAME='Dihedral RMSD test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
dihrms out dihrms.dat phi psi
EOF

  RunCpptraj "$UNITNAME"
  #DoTest dihedral.dat multidih.dat
fi

Disable() {
# Test nucleic acid chi
UNITNAME='Multidihedral Nucleotide CHI test'
CheckFor maxthreads 3
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm ../adh026.3.pdb 
trajin ../adh026.3.pdb 
multidihedral chin out chin.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest chin.dat.save chin.dat
fi

# Test protein dihedrals 
UNITNAME='Multidihedral protein dihedrals test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm ARG.mol2
trajin ARG.mol2
multidihedral out arg.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest arg.dat.save arg.dat
fi

UNITNAME='Multidihedral cyclic molecule test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm cyclic.mol2
trajin cyclic.mol2
multidihedral phi psi omega out cyclic.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest cyclic.dat.save cyclic.dat
fi
}
EndTest

exit 0
