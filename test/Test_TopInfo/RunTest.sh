#!/bin/bash

. ../MasterTest.sh

CleanFiles info.in atoms.dat residues.dat bonds.dat angles.dat dihedrals.dat

INPUT="-i info.in"
cat > info.in <<EOF
parm ../tz2.parm7
atoms :3
atoms :3 out atoms.dat

resinfo
resinfo out residues.dat
resinfo short
resinfo short out residues.dat

bonds :1
bonds :1 out bonds.dat
bonds @%N3 @%H
bonds @%N3 @%H out bonds.dat

angles @1
angles @1 out angles.dat
angles @%H @%N3 @%CT
angles @%H @%N3 @%CT out angles.dat

dihedralinfo @1
dihedralinfo @1 out dihedrals.dat
dihedralinfo @N @CA @CB @%H1
dihedralinfo @N @CA @CB @%H1 out dihedrals.dat
quit
EOF
RunCpptraj "Atom/bond info test."
DoTest atoms.dat.save atoms.dat
DoTest residues.dat.save residues.dat
DoTest bonds.dat.save bonds.dat
DoTest angles.dat.save angles.dat
DoTest dihedrals.dat.save dihedrals.dat

EndTest
exit 0
