#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in dssp.dat dssp.dat.sum dssp.sum.agr dssp.gnu

# Test 1
CheckNetcdf
cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc 
#secstruct out dssp.dat sumout dssp.sum.agr
secstruct out dssp.gnu sumout dssp.sum.agr 
EOF
INPUT="-i cpptraj.in"
TOP="../DPDP.parm7"
RunCpptraj "Secstruct (DSSP) command test."
DoTest dssp.gnu.save dssp.gnu
DoTest dssp.sum.agr.save dssp.sum.agr
CheckTest

# Test 2
CheckNetcdf
cat > cpptraj.in <<EOF
trajin ../DPDP.nc 
secstruct :10-22 out dssp.dat ptrajformat 
EOF
INPUT="-i cpptraj.in"
TOP="../DPDP.parm7"
RunCpptraj "Secstruct (DSSP) command test, Ptraj Format."
DoTest dssp.dat.save dssp.dat
DoTest dssp.dat.sum.save dssp.dat.sum
CheckTest


EndTest

exit 0
