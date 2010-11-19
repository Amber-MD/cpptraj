#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles goodrmsd.dat badrmsd.dat 

# Test 1
INPUT="-i goodtraj.in"
RunCpptraj "Broken Traj: Running good trajectory."

# Test 2
INPUT="-i brokentraj.in"
RunCpptraj "Broken Traj: Running broken trajectory."

DoTest goodrmsd.dat badrmsd.dat

CheckTest

EndTest

exit 0
