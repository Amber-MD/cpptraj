#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles goodrmsd.dat badrmsd.dat goodtraj.in brokentraj.in zip.gz zip.in ziprmsd.dat
RequiresNotParallel "Broken Traj"

# Test 1
cat > goodtraj.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.crd 1 100 
rmsd :2-18@N,CA,C out goodrmsd.dat
EOF
INPUT="-i goodtraj.in"
RunCpptraj "Broken Traj: Running good trajectory."

# Test 2
cat > brokentraj.in <<EOF
noprogress
parm ../tz2.parm7
trajin broken.tz2.crd
rmsd :2-18@N,CA,C out badrmsd.dat
EOF
INPUT="-i brokentraj.in"
RunCpptraj "Broken Traj: Running broken trajectory."
DoTest goodrmsd.dat badrmsd.dat

# Test 3
CheckZlib
cp broken.tz2.crd zip
gzip zip
cat > zip.in <<EOF
noprogress
parm ../tz2.parm7
trajin zip.gz
rmsd :2-18@N,CA,C out ziprmsd.dat
EOF
INPUT="-i zip.in"
RunCpptraj "Broken Traj: Running compressed broken trajectory."
DoTest goodrmsd.dat ziprmsd.dat

EndTest

exit 0
