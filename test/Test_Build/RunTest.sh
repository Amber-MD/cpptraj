#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in tz2.parm7 tz2.rst7

TESTNAME='Build tests'

INPUT='-i cpptraj.in'

if [ ! -z "$AMBERHOME" ] ; then
  echo "Warning: AMBERHOME is set; temporarily unsetting for $TESTNAME"
  unset AMBERHOME
fi

cat > cpptraj.in <<EOF
source leaprc.protein.ff14SB
build ../tz2.pdb parmout tz2.parm7 crdout tz2.rst7
EOF
RunCpptraj "$TESTNAME"
DoTest tz2.parm7.save tz2.parm7 -I %VERSION
DoTest tz2.rst7.save tz2.rst7

EndTest
