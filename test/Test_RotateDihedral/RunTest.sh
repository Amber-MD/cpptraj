#!/bin/bash

. ../MasterTest.sh

CleanFiles rotate.in tz2.rotate.1.mol2

INPUT='-i rotate.in'

# Rotate single dihedral to target value
cat > rotate.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc 1 1 name TZ2
rotatedihedral crdset TZ2 value 35 res 8 type chip
crdout TZ2 tz2.rotate.1.mol2
EOF
RunCpptraj "Rotate dihedral to target value."
DoTest tz2.rotate.1.mol2.save tz2.rotate.1.mol2

EndTest
exit 0
