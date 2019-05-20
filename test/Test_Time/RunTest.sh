#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in  withtime.rst7

INPUT='-i cpptraj.in'

TESTNAME='Set Time tests'
Requires maxthreads 1

UNITNAME='Set Time test'
cat > cpptraj.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
time time0 10.1
trajout withtime.rst7
go
EOF
RunCpptraj "$UNITNAME"
DoTest withtime.rst7.save withtime.rst7

EndTest
exit 0
