#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles phi psi omega dist_end_to_end.list test.mdcrd watershell.list
 
# Test 1
INPUT="ptraj.in"
TOP="compound.prmtop"
RunCpptraj "Comprehensive tests."
DoTest cpptraj.phi.save phi
DoTest cpptraj.psi.save psi
DoTest cpptraj.omega.save omega
DoTest cpptraj.dist_end_to_end.list.save dist_end_to_end.list
DoTest cpptraj.test.mdcrd.save test.mdcrd
DoTest watershell.list.save watershell.list

CheckTest

EndTest

exit 0
