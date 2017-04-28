#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles reordered.mol2 remap.in

INPUT="-i remap.in"
MaxThreads 1 "Remap Test"
if [ "$?" -eq 0 ] ; then
  # Test 1
  cat > remap.in <<EOF
readdata ../Test_AtomMap/atommap.dat.save name MyMap
parm ../Test_AtomMap/xtallig.mol2
trajin ../Test_AtomMap/xtallig.mol2
remap data MyMap:1
trajout reordered.mol2
EOF
  RunCpptraj "Remap Test"
  DoTest ../Test_AtomMap/reordered.mol2.save reordered.mol2
fi

EndTest

exit 0
