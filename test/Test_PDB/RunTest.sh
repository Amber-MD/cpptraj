#!/bin/bash

. ../MasterTest.sh

CleanFiles pdb.in test.pdb tz2.pqr.gb.pdb tz2.pqr.parse.pdb tz2.pqr.vdw.pdb

INPUT="-i pdb.in"

# Test read/write of residue numbers, insertion / altloc codes, etc
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

# Test writing PQR files with various radii options
NotParallel "PQR file write with various radii"
if [[ $? -eq 0 ]] ; then
  cat > pdb.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
trajout tz2.pqr.gb.pdb dumpq    # GB radii
trajout tz2.pqr.parse.pdb parse # PARSE radii
trajout tz2.pqr.vdw.pdb dumpr*  # VDW radii
EOF
  RunCpptraj "PQR file write with various radii"
  DoTest tz2.pqr.gb.pdb.save tz2.pqr.gb.pdb
  DoTest tz2.pqr.parse.pdb.save tz2.pqr.parse.pdb
  DoTest tz2.pqr.vdw.pdb.save tz2.pqr.vdw.pdb
fi


EndTest
exit 0
