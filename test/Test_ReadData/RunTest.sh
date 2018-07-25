#!/bin/bash

. ../MasterTest.sh

CleanFiles vector.in v6and7.dat rex-d.dat MD.ene.dat out.dx out.dat

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

UNITNAME='Read Amber output test'
cat > vector.in <<EOF
readdata md.initial.out md.restart.out name MD
writedata MD.ene.dat MD[*] prec 14.4
EOF
RunCpptraj "$UNITNAME"
DoTest MD.ene.dat.save MD.ene.dat

UNITNAME='Read CHARMM output test'
cat > vector.in <<EOF
readdata rex-d.out_0 name ENE
writedata rex-d.dat ENE[*]
EOF
RunCpptraj "$UNITNAME"
DoTest rex-d.dat.save rex-d.dat

# Make sure we can read data in and write it out properly
UNITNAME='3D data read/write test'
cat > vector.in <<EOF
readdata ../Test_Grid/out.dx.save name grid
writedata out.dx grid
EOF
RunCpptraj "$UNITNAME"
DoTest ../Test_Grid/out.dx.save out.dx

# Read data as DX, write standard
UNITNAME='Standard 3D data write test'
cat > vector.in <<EOF
readdata ../Test_Grid/out.dx.save name grid
writedata out.dat grid
EOF
RunCpptraj "$UNITNAME"

EndTest
  
exit 0
