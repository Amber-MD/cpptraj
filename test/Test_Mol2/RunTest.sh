#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles mol2.in L01.mol2 tz2.mol2 test.mol2 test1.mol2 test2.mol2

INPUT="-i mol2.in"

cat > mol2.in <<EOF
noprogress
parm charged.mol2
trajin charged.mol2 
trajout L01.mol2 mol2 title MOL

parm ../tz2.parm7
trajin ../tz2.rst7 parm tz2.parm7
trajout tz2.mol2 mol2 parm tz2.parm7
EOF
RunCpptraj "Mol2 Parm/Traj Read/Write Test."
DoTest L01.mol2.1.save L01.mol2
DoTest tz2.mol2.1.save tz2.mol2

cat > mol2.in <<EOF
parm ../tz2.pdb
trajin ../tz2.pdb
trajout test.mol2 mol2
EOF
RunCpptraj "PDB => Mol2"
DoTest test.mol2.save test.mol2

cat > mol2.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
trajout test1.mol2
EOF
RunCpptraj "Amber Top/Rst => Mol2"
DoTest test1.mol2.save test1.mol2

# SYBYL atom type conversion requires data in AMBERHOME
if [[ ! -z $AMBERHOME ]] ; then
  cat > mol2.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
trajout test2.mol2 sybyltype
EOF
  RunCpptraj "Amber Top/Rst => Mol2, SYBYL atom types"
  DoTest test2.mol2.save test2.mol2
else
  echo "Amber to SYBYL atom type conversion test requires AMBERHOME be set."
fi

EndTest

exit 0
