#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in dssp.dat dssp.dat.sum dssp.sum.agr dssp.gnu \
           dssp2.gnu dssp2.sum.agr total.agr assign.4.dat

CheckNetcdf
INPUT="-i cpptraj.in"
TOP="../DPDP.parm7"

# Test 1
cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc 
secstruct out dssp.gnu sumout dssp.sum.agr totalout total.agr 
EOF
RunCpptraj "Secstruct (DSSP) command test."
DoTest dssp.gnu.save dssp.gnu
DoTest dssp.sum.agr.save dssp.sum.agr
DoTest total.agr.save total.agr
CheckTest

# Test 2
cat > cpptraj.in <<EOF
trajin ../DPDP.nc 
secstruct :10-22 out dssp.dat ptrajformat 
EOF
RunCpptraj "Secstruct (DSSP) command test, Ptraj Format."
DoTest dssp.dat.save dssp.dat
DoTest dssp.dat.sum.save dssp.dat.sum
CheckTest

# Test 3
cat > cpptraj.in <<EOF
noprogress
parm ../tz2.nhe.parm7
trajin ../tz2.nhe.nc parmindex 1
trajin ../DPDP.nc
secstruct out dssp2.gnu sumout dssp2.sum.agr nostring
EOF
RunCpptraj "Secstruct (DSSP) command with changing number of residues."
DoTest dssp2.gnu.save dssp2.gnu
DoTest dssp2.sum.agr.save dssp2.sum.agr

# Test 4
TOP=""
cat > cpptraj.in <<EOF
parm test.4.pdb
trajin test.4.pdb
dssp N4 assignout assign.4.dat
EOF
RunCpptraj "Secstruct (DSSP): SS assign output test"
DoTest assign.4.dat.save assign.4.dat

EndTest

exit 0
