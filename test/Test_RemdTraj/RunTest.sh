#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles remd.in d1.offset.dat d1.crd.dat d1.nc.dat temp.crd.* \
           temp0.crd.* d1.ensemble.dat d1.ensemble.dat.? all.dat 

TESTNAME='Replica exchange trajectory tests'
Requires maxthreads 10

INPUT="-i remd.in"

# Test 0
UNITNAME='CRD Replica Trajectory Run with offset'
CheckFor maxthreads 5
if [ $? -eq 0 ] ; then
  cat > remd.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7 
trajin rem.crd.000 remdtraj remdtrajtemp 492.2 1 11 2
distance d1 out d1.offset.dat @1 @21
EOF
  RunCpptraj "$UNITNAME"
  DoTest d1.offset.dat.save d1.offset.dat
fi

# Test 1
cat > remd.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7 
trajin rem.crd.000 remdtraj remdtrajtemp 492.2 
distance d1 out d1.crd.dat @1 @21
EOF
RunCpptraj "CRD Replica Trajectory Run"
DoTest d1.crd.dat.save d1.crd.dat

# Test 2
UNITNAME='NetCDF Replica Trajectory Run test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > remd.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7 
trajin rem.nc.000 remdtraj remdtrajtemp 492.2
distance d1 out d1.nc.dat @1 @21
EOF
  RunCpptraj "$UNITNAME"
  DoTest d1.nc.dat.save d1.nc.dat
fi

# Remdout test
UNITNAME='CRD Replica Trajectory Run with remdout'
CheckFor nthreads 4
if [ $? -eq 0 ] ; then
  # Create trajectories at all temperatures.
  for T in 300.00 384.30 492.20 630.50 ; do
    cat > remd.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7 
trajin rem.crd.000 remdtraj remdtrajtemp $T
trajout temp0.crd.$T
EOF
  RunCpptraj "CRD Replica Trajectory Run: Generating $T traj"
  done
  # Convert ensemble to temperature trajectories in 1 step
  cat > remd.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7 
ensemble rem.crd.000
trajout temp.crd 
distance d1 out d1.ensemble.dat @1 @21
EOF
  RunCpptraj "CRD Replica Trajectory Run with remdout"
  if [ -z "$DO_PARALLEL" ] ; then
    DoTest d1.ensemble.dat.save d1.ensemble.dat
  else
    cat d1.ensemble.dat.? > all.dat
    DoTest all.dat.save all.dat
  fi
  DoTest temp0.crd.300.00 temp.crd.0
  DoTest temp0.crd.384.30 temp.crd.1
  DoTest temp0.crd.492.20 temp.crd.2
  DoTest temp0.crd.630.50 temp.crd.3
fi

EndTest

exit 0
