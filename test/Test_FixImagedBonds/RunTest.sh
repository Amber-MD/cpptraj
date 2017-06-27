#!/bin/bash

. ../MasterTest.sh

CleanFiles fix.in fixed.rst7 fixed.pdb unimage.crd

INPUT='-i fix.in'

MaxThreads 1 "Fix imaged bonds test (Ortho.)."
if [ $? -eq 0 ] ; then
  cat > fix.in <<EOF
parm MET.pdb
trajin MET.pdb
fiximagedbonds
trajout fixed.rst7
EOF
  RunCpptraj "Fix imaged bonds test."
  DoTest fixed.rst7.save fixed.rst7
fi

MaxThreads 10 "Fix imaged bonds test (Non-ortho.)."
if [ $? -eq 0 ] ; then
  cat > fix.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
strip :WAT
image byatom
fiximagedbonds
trajout unimage.crd
EOF
  RunCpptraj "Fix imaged bonds test (Non-ortho.)."
  DoTest unimage.crd.save unimage.crd
fi
EndTest
exit 0
