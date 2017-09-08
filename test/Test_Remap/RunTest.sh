#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles reordered.mol2 remap.in

TESTNAME='Remap test'
Requires maxthreads 1

INPUT="-i remap.in"
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

EndTest

exit 0
