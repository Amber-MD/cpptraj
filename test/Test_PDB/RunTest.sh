#!/bin/bash

. ../MasterTest.sh

CleanFiles pdb.in test.pdb tz2.pqr.gb.pdb tz2.pqr.parse.pdb \
           tz2.pqr.vdw.pdb chainA.dat oresnum.dat tz2.plain.pdb \
           2b5t.fromparm.pdb

INPUT="-i pdb.in"

TESTNAME='PDB format tests'
Requires maxthreads 1

# Test read/write of residue numbers, insertion / altloc codes, etc
Test1() {
  UNITNAME='PDB format read/write test'
cat > pdb.in <<EOF
parm 2b5t.pdb noconect
resinfo ::A out chainA.dat
resinfo :;2 out oresnum.dat
trajin 2b5t.pdb
trajout test.pdb teradvance sg "P 1"
EOF
  RunCpptraj "$UNITNAME"
  DoTest test.pdb.save test.pdb
  DoTest chainA.dat.save chainA.dat
  DoTest oresnum.dat.save oresnum.dat
}

# Test writing PQR files with various radii options
Test2() {
  UNITNAME='PQR file write with various radii'
  CheckFor notparallel
  if [ $? -eq 0 ] ; then
    cat > pdb.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
trajout tz2.pqr.gb.pdb dumpq    # GB radii
trajout tz2.pqr.parse.pdb parse # PARSE radii
trajout tz2.pqr.vdw.pdb dumpr*  # VDW radii
# Test default chain ID write
trajout tz2.plain.pdb pdbres
EOF
    RunCpptraj "$UNITNAME"
    DoTest tz2.pqr.gb.pdb.save tz2.pqr.gb.pdb
    DoTest tz2.pqr.parse.pdb.save tz2.pqr.parse.pdb
    DoTest tz2.pqr.vdw.pdb.save tz2.pqr.vdw.pdb
    DoTest tz2.plain.pdb.save tz2.plain.pdb
  fi
}

Test1
Test2

# Test adding back PDB info
UNITNAME='Adding PDB info to topology'
cat > pdb.in <<EOF
parm 2b5t.segment.pdb
trajin 2b5t.segment.pdb
change oresnums of :1 min 186 max 186
change oresnums of :6 min 187 max 187
change icodes of :2-5 min A max D resnum 186
change chainid of :1-6 to B
trajout 2b5t.fromparm.pdb
EOF
RunCpptraj "$UNITNAME"
DoTest 2b5t.fromparm.pdb.save 2b5t.fromparm.pdb

EndTest
exit 0
