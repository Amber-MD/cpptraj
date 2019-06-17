#!/bin/bash

. ../MasterTest.sh

CleanFiles dih.in dihrms.dat toref.dat previous.dat totraj.dat

TESTNAME='Dihedral RMSD tests'
Requires maxthreads 10

INPUT="-i dih.in"
UNITNAME='Dihedral RMSD test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
reference ../tz2.nc 1 [MyRef]
dihrms ToFirst out dihrms.dat noheader phi psi
dihrms ToRef ref [MyRef] out toref.dat noheader phi psi
dihrms ToPrev previous out previous.dat phi psi
# NOTE: All calculated values for same traj should be 0.0
dihrms ToTraj reftraj ../tz2.nc out totraj.dat phi psi
EOF
  RunCpptraj "$UNITNAME"
  DoTest dihrms.dat.save dihrms.dat
  DoTest dihrms.dat.save toref.dat
  DoTest previous.dat.save previous.dat
  DoTest totraj.dat.save totraj.dat
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
