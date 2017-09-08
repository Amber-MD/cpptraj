#!/bin/bash

. ../MasterTest.sh

CleanFiles fix.in fixed.rst7 fixed.pdb unimage.crd

INPUT='-i fix.in'

TESTNAME='Fix imaged bonds tests'
Requires maxthreads 10

UNITNAME='Fix imaged bonds test (Ortho.)'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > fix.in <<EOF
parm MET.pdb
trajin MET.pdb
fiximagedbonds
trajout fixed.rst7
EOF
  RunCpptraj "$UNITNAME"
  DoTest fixed.rst7.save fixed.rst7
fi

UNITNAME='Fix imaged bonds test (Non-ortho.)'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > fix.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
image byatom
fiximagedbonds :1-13
strip :WAT
trajout unimage.crd
EOF
  RunCpptraj "$UNITNAME"
  DoTest unimage.crd.save unimage.crd
fi
EndTest
exit 0
