#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles remd.in d1.offset.dat d1.crd.dat d1.nc.dat temp.crd.* temp0.crd.* d1.ensemble.dat d1.ensemble.dat.? all.dat 

if [[ -z $DO_PARALLEL ]] ; then
  CPPTRAJ_SERIAL=$CPPTRAJ
  CPPTRAJ_PARALLEL=$CPPTRAJ
  PARALLEL_CMD=""
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
CheckTest

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
CheckTest

# Test 1
CheckNetcdf
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
else
  echo "DO_PARALLEL is set. Skipping 'trajin' tests."
  CPPTRAJ_PARALLEL=$CPPTRAJ
  CPPTRAJ_SERIAL=${CPPTRAJ%".MPI"}
  if [[ ! -e $CPPTRAJ_SERIAL ]] ; then
    echo "Error: $CPPTRAJ_SERIAL not found."
    exit 1
  fi
  PARALLEL_CMD=$DO_PARALLEL
  DO_PARALLEL=""
fi

# Remdout test
# Create traj at all temperatures. Must be run in serial.
CPPTRAJ=$CPPTRAJ_SERIAL
INPUT="-i remd.in"
for T in 300.00 384.30 492.20 630.50 ; do
  cat > remd.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7 
trajin rem.crd.000 remdtraj remdtrajtemp $T
trajout temp0.crd.$T
EOF
RunCpptraj "CRD Replica Trajectory Run: Generating $T traj"
done

# Convert to temperature traj in 1 step
DO_PARALLEL=$PARALLEL_CMD
CPPTRAJ=$CPPTRAJ_PARALLEL
cat > remd.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7 
ensemble rem.crd.000 remdtraj remdtrajtemp 492.20 
trajout temp.crd 
distance d1 out d1.ensemble.dat @1 @21
EOF
INPUT="-i remd.in"
RunCpptraj "CRD Replica Trajectory Run with remdout"
if [[ -z $DO_PARALLEL ]] ; then
  DoTest d1.ensemble.dat.save d1.ensemble.dat
else
  cat d1.ensemble.dat.? > all.dat
  DoTest all.dat.save all.dat
fi
DoTest temp0.crd.300.00 temp.crd.0
DoTest temp0.crd.384.30 temp.crd.1
DoTest temp0.crd.492.20 temp.crd.2
DoTest temp0.crd.630.50 temp.crd.3
CheckTest

EndTest

exit 0
