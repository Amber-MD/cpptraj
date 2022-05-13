#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in sorted.remlog.crd.?

TESTNAME='H-REMD sorting tests'
Requires netcdf

INPUT='-i cpptraj.in'

UNITNAME='Sort by crdidx via remlog'
cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
ensemble rem.crd.001 remlog rem.log nstlim 1000 ntwx 1000
trajout sorted.remlog.crd
EOF
RunCpptraj "$UNITNAME"

EndTest
