#!/bin/bash

. ../MasterTest.sh

CleanFiles rotate.in tz2.rotate.?.mol2

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

# Rotate single dihedral by increment
cat > rotate.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc 1 1 name TZ2
rotatedihedral crdset TZ2 increment -19.3431 :8@N :8@CA :8@CB :8@CG
crdout TZ2 tz2.rotate.2.mol2
EOF
RunCpptraj "Rotate dihedral by increment"
DoTest tz2.rotate.1.mol2.save tz2.rotate.1.mol2

# Rotate as a TRAJ data set
cat > rotate.in <<EOF
parm ../tz2.parm7
loadtraj ../tz2.nc name TZ2
rotatedihedral crdset TZ2 value 35 :8@N :8@CA :8@CB :8@CG \
    name OUT
crdout OUT tz2.rotate.3.mol2
EOF
RunCpptraj "Rotate dihedral in TRAJ data set"
DoTest tz2.rotate.1.mol2.save tz2.rotate.3.mol2

EndTest
exit 0
