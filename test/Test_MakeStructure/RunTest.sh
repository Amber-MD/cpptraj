#!/bin/bash

. ../MasterTest.sh

CleanFiles ms.in pp2.rst7 hairpin.rst7 dihedrals.dat dihedrals?.dat fromref.rst7 fromref.pdb.1
INPUT="-i ms.in"

TESTNAME='Make structure tests'
Requires maxthreads 1

# Tests
MS1() {
  UNITNAME='Basic makestructure test'
  CheckFor netcdf
  if [ $? -eq 0 ] ; then
    cat > ms.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 1
# Tests SS arg
makestructure pp2:1-13
outtraj pp2.rst7 time0 1.0

# Tests SS arg and custom turn arg
makestructure extended:1,12 custom1:2-5:-80.0:130.0:-130.0:140.0 typeI':6-7 custom2:8-11:-140.0:170.0:-100.0:140.0 
multidihedral phi psi resrange 1-13 out dihedrals.dat

# Tests custom dih arg and custom ss arg
makestructure customdih:5:phi:90 custom:6,7:-70:60
multidihedral phi psi resrange 5-7 out dihedrals2.dat

# Test custom dih definition
makestructure chi1:8:N:CA:CB:CG:35
dihedral chi1 :8@N :8@CA :8@CB :8@CG out dihedrals3.dat
EOF
    RunCpptraj "$UNITNAME"
    DoTest pp2.rst7.save pp2.rst7
    DoTest dihedrals.dat.save dihedrals.dat -a 0.00001
    DoTest dihedrals2.dat.save dihedrals2.dat
    DoTest dihedrals3.dat.save dihedrals3.dat
  fi
}

Ref() {
# Test SS from reference
  cat > ms.in <<EOF
parm ../tz2.parm7
reference ../tz2.rst7
trajin pp2.rst7.save

makestructure "ref:1-13:tz2.rst7"
rmsd reference 
trajout fromref.pdb multi
EOF
  RunCpptraj "Makestructure test with reference."
  #DoTest ../tz2.rst7 fromref.rst7
  DoTest fromref.pdb.save fromref.pdb.1
}

MS1
Ref

EndTest
exit 0
