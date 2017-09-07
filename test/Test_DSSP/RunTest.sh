#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in dssp.dat dssp.dat.sum dssp.sum.agr dssp.gnu \
           dssp2.gnu dssp2.sum.agr total.agr assign.4.dat \
           DSSP.assign.dat

INPUT="-i cpptraj.in"
TOP="../DPDP.parm7"

# Test 1
CheckNetcdf "DSSP tests"
if [ $? -ne 0 ] ; then
  SkipCheck "DSSP tests"
else
  cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc 
secstruct out dssp.gnu sumout dssp.sum.agr totalout total.agr \
          assignout DSSP.assign.dat
EOF
  RunCpptraj "Secstruct (DSSP) command test."
  DoTest dssp.gnu.save dssp.gnu
  DoTest dssp.sum.agr.save dssp.sum.agr
  DoTest total.agr.save total.agr
  DoTest DSSP.assign.dat.save DSSP.assign.dat

  # Test 2
  cat > cpptraj.in <<EOF
trajin ../DPDP.nc 
secstruct :10-22 out dssp.dat ptrajformat 
EOF
  RunCpptraj "Secstruct (DSSP) command test, Ptraj Format."
  DoTest dssp.dat.save dssp.dat
  DoTest dssp.dat.sum.save dssp.dat.sum

  # Test 3
  TITLE="Secstruct (DSSP) command with changing number of residues."
  NotParallel "$TITLE"
  if [ $? -ne 0 ] ; then
    SkipCheck "$TITLE"
  else
    cat > cpptraj.in <<EOF
noprogress
parm ../tz2.nhe.parm7
trajin ../tz2.nhe.nc parmindex 1
trajin ../DPDP.nc
secstruct out dssp2.gnu sumout dssp2.sum.agr nostring
EOF
    RunCpptraj "$TITLE"
    DoTest dssp2.gnu.save dssp2.gnu
    DoTest dssp2.sum.agr.save dssp2.sum.agr
  fi
fi # END if netcdf

# Test 4
TESTNAME='Secstruct (DSSP): SS assign output test'
MaxThreads 1 "$TESTNAME"
if [ $? -ne 0 ] ; then
  SkipCheck "$TESTNAME"
else
  TOP=""
  cat > cpptraj.in <<EOF
parm test.4.pdb
trajin test.4.pdb
dssp N4 assignout assign.4.dat
EOF
  RunCpptraj "Secstruct (DSSP): SS assign output test"
  DoTest assign.4.dat.save assign.4.dat
fi

EndTest

exit 0
