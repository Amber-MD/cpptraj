#!/bin/bash

. ../MasterTest.sh

CleanFiles scale.in scale.res1.0.5.dat scale.N.C.CA.0.5.dat

INPUT='-i scale.in'

cat > scale.in <<EOF
parm ../tz2.parm7
scaledihedralk 0.5
dihedrals :1 out scale.res1.0.5.dat
EOF
RunCpptraj "Scale dihedral force constant test"
DoTest scale.res1.0.5.dat.save scale.res1.0.5.dat

cat > scale.in <<EOF
parm ../tz2.parm7
scaledihedralk 0.5 @N,C,CA useall
dihedrals :1 out scale.N.C.CA.0.5.dat
EOF
RunCpptraj "Scale dihedral force constant with 'useall'"
DoTest scale.N.C.CA.0.5.dat.save scale.N.C.CA.0.5.dat

EndTest
exit 0
