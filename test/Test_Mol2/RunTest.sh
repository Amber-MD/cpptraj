#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles mol2.in L01.mol2.1 trpcage.mol2.1 

# Test 1
cat > mol2.in <<EOF
noprogress
parm charged.mol2
trajin charged.mol2 
trajout L01.mol2 mol2 title MOL

parm ../trpcage.parm7
trajin ../trpcage.rst7 parm trpcage.parm7
trajout trpcage.mol2 mol2 parm trpcage.parm7
EOF
INPUT="-i mol2.in"
RunCpptraj "Mol2 Parm/Traj Read/Write Test."
DoTest L01.mol2.1.save L01.mol2.1
DoTest trpcage.mol2.1.save trpcage.mol2.1
CheckTest

EndTest

exit 0
