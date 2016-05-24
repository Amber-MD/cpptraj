#!/bin/bash

. ../MasterTest.sh

CleanFiles pdb.in test.pdb

INPUT="-i pdb.in"

MaxThreads 1 "PDB format read/write test."
if [[ $? -eq 0 ]] ; then
  cat >> pdb.in <<EOF
parm 2b5t.pdb noconect
trajin 2b5t.pdb
trajout test.pdb teradvance sg "P 1"
EOF
  RunCpptraj "PDB format read/write test."
  DoTest test.pdb.save test.pdb
fi

EndTest
exit 0
