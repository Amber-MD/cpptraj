#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles remd.in d1.offset.dat d1.crd.dat d1.nc.dat

# Test 0
cat > remd.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7 
trajin rem.crd.000 remdtraj remdtrajtemp 492.2 1 11 2
distance d1 out d1.offset.dat @1 @21
EOF
INPUT="-i remd.in"
RunCpptraj "CRD Replica Trajectory Run with offset"
DoTest d1.offset.dat.save d1.offset.dat
#CheckTest

# Test
cat > remd.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7 
trajin rem.crd.000 remdtraj remdtrajtemp 492.2 
distance d1 out d1.crd.dat @1 @21
EOF
INPUT="-i remd.in"
RunCpptraj "CRD Replica Trajectory Run"
DoTest d1.crd.dat.save d1.crd.dat

# Test 1
cat > remd.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7 
trajin rem.nc.000 remdtraj remdtrajtemp 492.2
distance d1 out d1.nc.dat @1 @21
EOF
INPUT="-i remd.in"
RunCpptraj "NETCDF Replica Trajectory Run"
DoTest d1.nc.dat.save d1.nc.dat
CheckTest

EndTest

exit 0
