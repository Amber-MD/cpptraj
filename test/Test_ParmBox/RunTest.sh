#!/bin/bash

. ../MasterTest.sh

CleanFiles parmbox.in out.parm7

INPUT='-i parmbox.in'

cat > parmbox.in <<EOF
parm ../Test_SymmRmsd/TYR.parm7
parmbox truncoct x 10.0
parmwrite out out.parm7
EOF
RunCpptraj "Parmbox test"
DoTest out.parm7.save out.parm7 -I %VERSION

EndTest
exit 0
