#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in unwrap.crd unwrap.ortho.crd
TESTNAME='Unwrap tests'
Requires netcdf notparallel

INPUT="ptraj.in"
TOP="../tz2.truncoct.parm7"
cat > ptraj.in <<EOF
trajin ../tz2.truncoct.nc 1 2
unwrap 
trajout unwrap.crd title "Test"
EOF
RunCpptraj "Unwrap non-orthogonal test"
DoTest unwrap.crd.save unwrap.crd

TOP="../tz2.ortho.parm7"
cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc 1 2
unwrap 
trajout unwrap.ortho.crd title "Test"
EOF
RunCpptraj "Unwrap orthogonal test"
DoTest unwrap.ortho.crd.save unwrap.ortho.crd

EndTest

exit 0
