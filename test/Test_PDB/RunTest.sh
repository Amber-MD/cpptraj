#!/bin/bash

. ../MasterTest.sh

CleanFiles pdb.in test.pdb tz2.pqr.gb.pdb tz2.pqr.parse.pdb tz2.pqr.vdw.pdb

INPUT="-i pdb.in"

TESTNAME='PDB format tests'
Requires maxthreads 1

# Test read/write of residue numbers, insertion / altloc codes, etc
UNITNAME='PDB format read/write test'
cat >> pdb.in <<EOF
parm 2b5t.pdb noconect
trajin 2b5t.pdb
trajout test.pdb teradvance sg "P 1"
EOF
RunCpptraj "$UNITNAME"
DoTest test.pdb.save test.pdb

# Test writing PQR files with various radii options
UNITNAME='PQR file write with various radii'
CheckFor notparallel
if [ $? -eq 0 ] ; then
  cat > pdb.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
trajout tz2.pqr.gb.pdb dumpq    # GB radii
trajout tz2.pqr.parse.pdb parse # PARSE radii
trajout tz2.pqr.vdw.pdb dumpr*  # VDW radii
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2.pqr.gb.pdb.save tz2.pqr.gb.pdb
  DoTest tz2.pqr.parse.pdb.save tz2.pqr.parse.pdb
  DoTest tz2.pqr.vdw.pdb.save tz2.pqr.vdw.pdb
fi

EndTest
exit 0
