#!/bin/bash

. ../MasterTest.sh

CleanFiles info.in atoms.dat

INPUT="-i info.in"
cat > info.in <<EOF
parm ../tz2.parm7
atoms :3
atoms :3 out atoms.dat
quit
EOF
RunCpptraj "Atom info test."
DoTest atoms.dat.save atoms.dat

EndTest
exit 0
