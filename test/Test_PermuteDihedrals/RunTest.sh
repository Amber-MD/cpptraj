#!/bin/bash

. ../MasterTest.sh

CleanFiles ds.in rotations.nc rotations.mdcrd random.mol2.1
INPUT="ds.in"
TOP=../tz2.parm7

MaxThreads 1 "PermuteDihedrals tests"
if [[ $? -ne 0 ]] ; then
  EndTest
  exit 0
fi

cat > ds.in <<EOF
reference ../tz2.rst7 [TZ2]
permutedihedrals crdset [TZ2] interval -120 outtraj rotations.mdcrd phi psi
EOF
RunCpptraj "PermuteDihedrals interval -120, phi psi"
#DoTest rotations.nc.save rotations.nc
DoTest rotations.mdcrd.save rotations.mdcrd

cat > ds.in <<EOF
reference ../tz2.rst7 [TZ2]
permutedihedrals crdset [TZ2] random rseed 1 check maxfactor 10 phi psi
crdout [TZ2] random.mol2 multi
EOF
RunCpptraj "PermuteDihedrals, random rotations"
DoTest random.mol2.save random.mol2.1

EndTest
exit 0

