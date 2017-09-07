#!/bin/bash

. ../MasterTest.sh

CleanFiles fix.in fixed.rst7 fixed.pdb unimage.crd

INPUT='-i fix.in'

RequiresMaxThreads 10 "Fix imaged bonds tests"

TESTNAME='Fix imaged bonds test (Ortho.)'
MaxThreads 1 "$TESTNAME"
if [ $? -ne 0 ] ; then
  SkipCheck "$TESTNAME"
else
  cat > fix.in <<EOF
parm MET.pdb
trajin MET.pdb
fiximagedbonds
trajout fixed.rst7
EOF
  RunCpptraj "$TESTNAME"
  DoTest fixed.rst7.save fixed.rst7
fi

TESTNAME='Fix imaged bonds test (Non-ortho.)'
CheckNetcdf "$TESTNAME"
if [ $? -ne 0 ] ; then
  SkipCheck "$TESTNAME"
else
  cat > fix.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
image byatom
fiximagedbonds :1-13
strip :WAT
trajout unimage.crd
EOF
  RunCpptraj "$TESTNAME"
  DoTest unimage.crd.save unimage.crd
fi
EndTest
exit 0
