#!/bin/bash

. ../MasterTest.sh

CleanFiles replicate.in cell.mol2 

TESTNAME='Replicate cell test'
Requires netcdf

INPUT='-i replicate.in'

cat > replicate.in <<EOF
parm ../Test_SymmRmsd/TYR.parm7
trajin ../Test_SymmRmsd/TYR.nc 1 1
box x 20 y 20 z 20 beta 90
replicatecell out cell.mol2 dir 100 dir 0-10 dir 001
EOF
RunCpptraj "$TESTNAME"
DoTest cell.mol2.save cell.mol2

EndTest
exit 0
