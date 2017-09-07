#!/bin/bash

. ../MasterTest.sh

CleanFiles compare.in parm.diff

INPUT="-i compare.in"
TESTNAME="Compare topology test"
RequiresNotParallel "$TESTNAME"

cat > compare.in <<EOF
parm Protein.ff14SB.parm7
parm Protein.ff99.parm7
comparetop parmindex 0 parmindex 1 out parm.diff
EOF
RunCpptraj "$TESTNAME"
DoTest parm.diff.save parm.diff

EndTest
exit 0
