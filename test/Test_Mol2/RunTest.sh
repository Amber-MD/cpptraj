#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles mol2.in L01.mol2 tz2.mol2 

# Test 1
cat > mol2.in <<EOF
noprogress
parm charged.mol2
trajin charged.mol2 
trajout L01.mol2 mol2 title MOL

parm ../tz2.parm7
trajin ../tz2.rst7 parm tz2.parm7
trajout tz2.mol2 mol2 parm tz2.parm7
EOF
INPUT="-i mol2.in"
RunCpptraj "Mol2 Parm/Traj Read/Write Test."
DoTest L01.mol2.1.save L01.mol2
DoTest tz2.mol2.1.save tz2.mol2
CheckTest

EndTest

exit 0
