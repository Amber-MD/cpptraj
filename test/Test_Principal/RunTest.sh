#!/bin/bash

# Test of principal 
. ../MasterTest.sh

CleanFiles principal.in Ctest.pdb principal.dat Ctest.crd

TOP="../Test_IRED/1IEE_A_prot.prmtop"
INPUT="principal.in"

cat > principal.in <<EOF
noprogress
trajin ../Test_IRED/1IEE_A_test.mdcrd 1 10
principal * dorotation mass out principal.dat
trajout Ctest.crd
EOF
RunCpptraj "Principal Test"
DoTest Ctest.crd.save Ctest.crd
DoTest principal.dat.save principal.dat

EndTest

exit 0
