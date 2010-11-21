#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles phi psi omega dist_end_to_end.list test.mdcrd trajectory.netcdf trajectory_test.mdcrd
 
# Test 1
INPUT="ptraj.in"
TOP="compound.prmtop"
RunCpptraj "PTRAJ comprehensive tests."
DoTest cpptraj.phi.save phi
DoTest cpptraj.psi.save psi
DoTest cpptraj.omega.save omega
DoTest cpptraj.dist_end_to_end.list.save dist_end_to_end.list
DoTest cpptraj.test.mdcrd.save test.mdcrd

CheckTest

# Test 2
CheckNetcdf
INPUT="ptraj_netcdf.in"
TOP="compound.prmtop"
RunCpptraj "PTRAJ comprehensive mdcrd -> netcdf test."
INPUT="ptraj_mdcrd.in"
RunCpptraj "PTRAJ comprehensive netcdf -> mdcrd test."
#DoTest trajectory.mdcrd trajectory_test.mdcrd
DoTest cpptraj.trajectory.mdcrd trajectory_test.mdcrd 
CheckTest

EndTest

exit 0
