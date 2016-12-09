#!/bin/bash

. ../MasterTest.sh

CleanFiles info.in atoms.dat bonds.dat

INPUT="-i info.in"
cat > info.in <<EOF
parm ../tz2.parm7
atoms :3
atoms :3 out atoms.dat

bonds :1
bonds :1 out bonds.dat
bonds @%N3 @%H
bonds @%N3 @%H out bonds.dat
quit
EOF
RunCpptraj "Atom/bond info test."
DoTest atoms.dat.save atoms.dat
DoTest bonds.dat.save bonds.dat

EndTest
exit 0
