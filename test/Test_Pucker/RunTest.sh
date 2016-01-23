#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles pucker.in pucker.dat Ptest.in CremerF.dat CremerP.dat

# Test 1
Nucleic() {
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
}

Furanoid() {
  TOP=""
  INPUT="-i Ptest.in"
  cat > Ptest.in <<EOF
parm Furanoid.mol2
trajin Furanoid.mol2
pucker Furanoid @C2 @C3 @C4 @C5 @O2 cremer out CremerF.dat amplitude range360
EOF
  RunCpptraj "5-member ring pucker, Cremer & Pople Furanoid test."
  DoTest CremerF.dat.save CremerF.dat
}

Pyranoid() {
  TOP=""
  INPUT="-i Ptest.in"
  cat > Ptest.in <<EOF
parm Pyranoid.mol2
trajin Pyranoid.mol2
pucker Pyranoid @C1 @C2 @C3 @C4 @C5 @O5 cremer out CremerP.dat amplitude theta range360
EOF
  RunCpptraj "6-member ring pucker, Cremer & Pople Pyranoid test."
  DoTest CremerP.dat.save CremerP.dat
}

MaxThreads 3 "Pucker command test"
if [[ $? -eq 0 ]] ; then
  Nucleic
fi
MaxThreads 1 "Pyranoid/furanoid pucker tests."
if [[ $? -eq 0 ]] ; then
  Furanoid
  Pyranoid
fi

EndTest

exit 0
