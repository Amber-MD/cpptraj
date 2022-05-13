#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in sorted.remlog.crd.? sorted.crdidx.crd.?

TESTNAME='H-REMD sorting tests'
Requires netcdf nthreads 4

INPUT='-i cpptraj.in'

UNITNAME='Sort by crdidx via remlog'
cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
ensemblesize 4
ensemble rem.crd.001 remlog rem.log nstlim 1000 ntwx 1000
trajout sorted.remlog.crd
EOF
RunCpptraj "$UNITNAME"
DoTest sorted.remlog.crd.0.save sorted.remlog.crd.0
DoTest sorted.remlog.crd.1.save sorted.remlog.crd.1
DoTest sorted.remlog.crd.2.save sorted.remlog.crd.2
DoTest sorted.remlog.crd.3.save sorted.remlog.crd.3

UNITNAME='Sort by crdidx from trajectory'
cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
ensemblesize 4
ensemble rem.crd.001 bycrdidx 
trajout sorted.crdidx.crd
EOF
RunCpptraj "$UNITNAME"
DoTest sorted.remlog.crd.0.save sorted.crdidx.crd.0
DoTest sorted.remlog.crd.1.save sorted.crdidx.crd.1
DoTest sorted.remlog.crd.2.save sorted.crdidx.crd.2
DoTest sorted.remlog.crd.3.save sorted.crdidx.crd.3

EndTest
