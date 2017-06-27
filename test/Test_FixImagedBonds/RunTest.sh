#!/bin/bash

. ../MasterTest.sh

CleanFiles fix.in fixed.rst7 fixed.pdb

INPUT='-i fix.in'

MaxThreads 1 "Fix imaged bonds test."
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
EndTest
exit 0
