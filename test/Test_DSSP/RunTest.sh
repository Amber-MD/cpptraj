#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in dssp.dat dssp.sum.agr

# Test 1
cat > cpptraj.in <<EOF
noprogress
trajin ../DPDP.nc 
secstruct out dssp.dat sumout dssp.sum.agr
EOF
INPUT="-i cpptraj.in"
TOP="../DPDP.mod.GA12.parm7"
RunCpptraj "Secstruct (DSSP) command test."
DoTest dssp.dat.save dssp.dat
DoTest dssp.sum.agr.save dssp.sum.agr
CheckTest

EndTest

exit 0
