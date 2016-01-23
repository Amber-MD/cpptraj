#!/bin/bash

# Test of principal 
. ../MasterTest.sh

CleanFiles principal.in Ctest.pdb principal.dat Ctest.crd eigen.dat

TOP="../Test_IRED/1IEE_A_prot.prmtop"
INPUT="principal.in"

if [[ -z $DO_PARALLEL ]] ; then
  cat > principal.in <<EOF
noprogress
trajin ../Test_IRED/1IEE_A_test.mdcrd 1 10
principal * dorotation mass out principal.dat name All
trajout Ctest.crd
run
writedata eigen.dat All[eval] All[evec]
EOF
  RunCpptraj "Principal Test"
  DoTest Ctest.crd.save Ctest.crd
  DoTest principal.dat.save principal.dat
  DoTest eigen.dat.save eigen.dat
else
  # 'out' not supported in parallel
  cat > principal.in <<EOF
noprogress
trajin ../Test_IRED/1IEE_A_test.mdcrd 1 10
principal * dorotation mass name All
trajout Ctest.crd
run
writedata eigen.dat All[eval] All[evec]
EOF
  RunCpptraj "Principal Test"
  DoTest Ctest.crd.save Ctest.crd
  DoTest eigen.dat.save eigen.dat
fi

EndTest

exit 0
