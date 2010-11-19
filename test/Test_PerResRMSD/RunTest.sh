#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles PerResRMSD.dat

# Test 1
INPUT="-i cpptraj.in"
RunCpptraj "Per-Residue RMSD Test."
DoTest PerResRMSD.dat.save PerResRMSD.dat
CheckTest

EndTest

exit 0
