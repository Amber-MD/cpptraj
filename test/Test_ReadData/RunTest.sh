#!/bin/bash

. ../MasterTest.sh

CleanFiles vector.in v6and7.dat rex-d.dat

TESTNAME='Read data tests'

INPUT="-i vector.in"
# Test read/append of vector dataset
UNITNAME='Read vector data test'
cat > vector.in <<EOF
readdata ../Test_Vector/vtest.dat.6.save vector name v6and7
readdata ../Test_Vector/vtest.dat.7.save vector name v6and7
writedata v6and7.dat v6and7
EOF
RunCpptraj "$UNITNAME"
DoTest v6and7.dat.save v6and7.dat

UNITNAME='Read CHARMM output test'
cat > vector.in <<EOF
readdata rex-d.out_0 name ENE
writedata rex-d.dat ENE[*]
EOF
RunCpptraj "$UNITNAME"
DoTest rex-d.dat.save rex-d.dat

EndTest
  
exit 0
