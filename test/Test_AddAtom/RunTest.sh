#!/bin/bash

. ../MasterTest.sh

CleanFiles addatom.in tz2.addatom.pdb tz2.mask.dat tz2.mg.mol2

TESTNAME='Add atom tests'

INPUT='-i addatom.in'

Requires maxthreads 1

cat > addatom.in <<EOF
parm ../tz2.pdb
trajin ../tz2.pdb 1 1
addatom aname TEMP rname TMP elt H xyz 1 1 1
mask :TMP<@3.0&!:TMP out tz2.mask.dat name TZ2
trajout tz2.addatom.pdb
run
EOF
RunCpptraj "Add atom test"
DoTest tz2.addatom.pdb.save tz2.addatom.pdb
DoTest tz2.mask.dat.save tz2.mask.dat

cat > addatom.in <<EOF
parm ../tz2.parm7
trajin ../tz2.pdb 1 1
strip !:1
addatom aname MG rname MG mass 24.3 charge 2 xyz 3 4 5
trajout tz2.mg.mol2
EOF
RunCpptraj "Add atom test, MG"
DoTest tz2.mg.mol2.save tz2.mg.mol2

EndTest
