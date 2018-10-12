#!/bin/bash

. ../MasterTest.sh

TESTNAME='Charmm parameters tests'

CleanFiles cpptraj.in test3.parm7

INPUT='-i cpptraj.in'

UNITNAME='Test updating topology with Charmm parameters (single residue)'
cat > cpptraj.in <<EOF
parm lys.parm7
readdata kcx.str as charmmrtfprm name MyParm
list dataset
updateparameters parmindex 0 setname MyParm
parmwrite out test3.parm7
EOF
RunCpptraj "$UNITNAME"
DoTest test3.parm7.save test3.parm7 -I %VERSION

EndTest
exit 0
