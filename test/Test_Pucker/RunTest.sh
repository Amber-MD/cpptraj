#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles pucker.in pucker.dat

# Test 1
TOP=../adh026.3.pdb
INPUT=pucker.in
cat > pucker.in <<EOF
trajin ../adh026.3.pdb
pucker p1-as :1@C1' :1@C2' :1@C3' :1@C4' :1@O4' out pucker.dat 
pucker p2-as :2@C1' :2@C2' :2@C3' :2@C4' :2@O4' out pucker.dat 
pucker p3-as :3@C1' :3@C2' :3@C3' :3@C4' :3@O4' out pucker.dat 
pucker p1-cp :1@C1' :1@C2' :1@C3' :1@C4' :1@O4' out pucker.dat cremer
pucker p2-cp :2@C1' :2@C2' :2@C3' :2@C4' :2@O4' out pucker.dat cremer
pucker p3-cp :3@C1' :3@C2' :3@C3' :3@C4' :3@O4' out pucker.dat cremer
EOF
RunCpptraj "Pucker command test"
DoTest pucker.dat.save pucker.dat

CheckTest

EndTest

exit 0
