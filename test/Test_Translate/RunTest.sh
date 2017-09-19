#!/bin/bash

. ../MasterTest.sh

CleanFiles translate.in translate.2.11.mol2

INPUT='-i translate.in'

cat > translate.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
strip !(:2,11)
translate :2 x  1.0 y  10.0 z -0.5
translate :1 x  2.0 y -10.0 z -0.5
trajout translate.2.11.mol2
EOF
RunCpptraj "Translate test"
DoTest translate.2.11.mol2.save translate.2.11.mol2

EndTest
exit 0
