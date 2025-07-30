#!/bin/bash

. ../MasterTest.sh

CleanFiles cif.in 1LE1.pdb temp?.crd Cpptraj.pdb

TESTNAME='CIF tests'
Requires maxthreads 6

INPUT="-i cif.in"
UNITNAME='CIF read test'
CheckFor notparallel
if [ $? -eq 0 ] ; then
  cat > cif.in <<EOF
parm 1LE1.cif
trajin 1LE1.cif
trajout 1LE1.pdb
EOF
  RunCpptraj "$UNITNAME"
  DoTest 1LE1.pdb.save 1LE1.pdb
fi

# Test reading in _chem_comp_atom block
UNITNAME='Read in _chem_comp_atom block'
CheckFor notparallel
if [ $? -eq 0 ] ; then
  cat > cif.in <<EOF
parm CRO.cif
trajin CRO.cif
trajout Cpptraj.pdb
EOF
  RunCpptraj "Read in _chem_comp_atom block"
  DoTest Cpptraj.pdb.save Cpptraj.pdb
fi

# Test read with offset
cat > cif.in <<EOF
parm 1LE1.pdb.save
trajin 1LE1.pdb.save 3 last 3
trajout temp1.crd
EOF
RunCpptraj "Generate offset traj from PDB"
cat > cif.in <<EOF
parm 1LE1.cif
trajin 1LE1.cif 3 last 3
trajout temp2.crd
EOF
RunCpptraj "Generate offset traj from CIF"
DoTest temp1.crd temp2.crd

EndTest
exit 0
