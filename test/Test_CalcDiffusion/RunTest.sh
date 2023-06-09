#!/bin/bash

. ../MasterTest.sh

CleanFiles diff.in unwrap.out unwrap.dcd tz2.diff.dat tz2.traj.diff.dat

INPUT='-i diff.in'

TESTNAME='Diffusion analysis tests'
Requires netcdf

UNITNAME='Diffusion analysis calculation test'
if [ -z "$DO_PARALLEL" ] ; then
  CheckFor maxthreads 1
  if [ $? -eq 0 ] ; then
    cat > diff.in <<EOF
parm ../tz2.ortho.parm7
loadcrd ../tz2.ortho.nc name TZ2 1 10 
crdaction TZ2 unwrap bymol
runanalysis calcdiffusion crdset TZ2 out tz2.diff.dat :WAT@O
EOF
    RunCpptraj "$UNITNAME"
    DoTest tz2.diff.dat.save tz2.diff.dat
  fi
else
  CheckFor maxthreads 6
  if [ $? -eq 0 ] ; then
    # First pass, unwrap in serial
    cat > diff.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
unwrap bymol
trajout unwrap.dcd
EOF
    $CPPTRAJ -i diff.in -o unwrap.out
    # Second pass, calc diffusion
    cat > diff.in <<EOF
parm ../tz2.ortho.parm7
loadcrd unwrap.dcd name TZ2 1 10
runanalysis calcdiffusion crdset TZ2 out tz2.diff.dat :WAT@O
EOF
    RunCpptraj "$UNITNAME"
    DoTest tz2.diff.dat.save tz2.diff.dat
    # Test with TRAJ coords
    cat > diff.in <<EOF
parm ../tz2.ortho.parm7
loadtraj unwrap.dcd name TZ2 1 10
runanalysis calcdiffusion crdset TZ2 out tz2.traj.diff.dat :WAT@O
EOF
    RunCpptraj "$UNITNAME, TRAJ set"
    DoTest tz2.diff.dat.save tz2.traj.diff.dat
  fi
fi

EndTest
exit 0
