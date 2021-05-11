#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles pucker.in nucleic.dat furanoid.dat pyranoid.type.dat pyranoid.auto.dat

TESTNAME='MultiPucker tests'
Requires maxthreads 3

# Test 1
Nucleic() {
  TOP=../adh026.3.pdb
  INPUT=pucker.in
  cat > pucker.in <<EOF
trajin ../adh026.3.pdb
multipucker ADHas resrange 1-3 altona out nucleic.dat
multipucker ADHcp resrange 1-3 cremer out nucleic.dat
#pucker p1-as :1@C1' :1@C2' :1@C3' :1@C4' :1@O4' out pucker.dat 
#pucker p2-as :2@C1' :2@C2' :2@C3' :2@C4' :2@O4' out pucker.dat 
#pucker p3-as :3@C1' :3@C2' :3@C3' :3@C4' :3@O4' out pucker.dat 
#pucker p1-cp :1@C1' :1@C2' :1@C3' :1@C4' :1@O4' out pucker.dat cremer
#pucker p2-cp :2@C1' :2@C2' :2@C3' :2@C4' :2@O4' out pucker.dat cremer
#pucker p3-cp :3@C1' :3@C2' :3@C3' :3@C4' :3@O4' out pucker.dat cremer
EOF
  RunCpptraj "MultiPucker command test"
  DoTest nucleic.dat.save nucleic.dat
}

Furanoid() {
  UNITNAME='MultiPucker 5-member ring pucker, Cremer & Pople Furanoid test'
  CheckFor maxthreads 1
  if [ $? -eq 0 ] ; then
    TOP=""
    INPUT="-i Ptest.in"
    cat > Ptest.in <<EOF
parm ../Test_Pucker/Furanoid.mol2
trajin ../Test_Pucker/Furanoid.mol2
multipucker Furanoid puckertype furanoid:C2:C3:C4:C5:O2 cremer \
  out furanoid.dat amplitude ampout furanoid.dat range360
EOF
    RunCpptraj "$UNITNAME"
    DoTest furanoid.dat.save furanoid.dat
  fi
}

Pyranoid() {
  UNITNAME='MultiPucker 6-member ring pucker, Cremer & Pople Pyranoid test'
  CheckFor maxthreads 1
  if [ $? -eq 0 ] ; then
    TOP=""
    INPUT="-i Ptest.in"
    cat > Ptest.in <<EOF
parm ../Test_Pucker/Pyranoid.mol2
trajin ../Test_Pucker/Pyranoid.mol2
multipucker Pyranoid puckertype pyranoid:C1:C2:C3:C4:C5:O5 cremer \
  out pyranoid.type.dat \
  amplitude ampout pyranoid.type.dat\
  theta thetaout pyranoid.type.dat range360
multipucker Pyr pyranose cremer \
  out pyranoid.auto.dat \
  amplitude ampout pyranoid.auto.dat \
  theta thetaout pyranoid.auto.dat range360
EOF
    RunCpptraj "$UNITNAME"
    DoTest pyranoid.type.dat.save pyranoid.type.dat
    DoTest pyranoid.auto.dat.save pyranoid.auto.dat
  fi
}

Nucleic
Furanoid
Pyranoid

EndTest

exit 0
