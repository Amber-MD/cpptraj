#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in tz2.parm7 tz2.rst7 tz2.coords.parm7 tz2.coords.rst7

TESTNAME='Build/Source tests'

INPUT='-i cpptraj.in'

if [ ! -z "$AMBERHOME" ] ; then
  echo "Warning: AMBERHOME is set; temporarily unsetting for $TESTNAME"
  unset AMBERHOME
fi

UNITNAME='Source and Build test, direct'
cat > cpptraj.in <<EOF
source leaprc.protein.ff14SB
build ../tz2.pdb parmout tz2.parm7 crdout tz2.rst7
EOF
RunCpptraj "$UNITNAME"
DoTest tz2.parm7.save tz2.parm7 -I %VERSION
DoTest tz2.rst7.save tz2.rst7

UNITNAME='Source and Build test, COORDS'
cat > cpptraj.in <<EOF
source leaprc.protein.ff14SB
readdata ../tz2.pdb as coords name MyCrd
build crdset MyCrd parmout tz2.coords.parm7 crdout tz2.coords.rst7
EOF
RunCpptraj "$UNITNAME"
DoTest tz2.parm7.save tz2.coords.parm7 -I %VERSION
DoTest tz2.rst7.save tz2.coords.rst7

EndTest
