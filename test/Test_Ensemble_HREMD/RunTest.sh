#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in sorted.remlog.crd.? sorted.crdidx.crd.? \
           distances.dat.? avg.distances.dat.? \
           distances.dat avg.distances.dat

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

# FIXME distance currently only sets up first ensemble member initially.
#       Should do all members.
UNITNAME='Ensemble analysis test'
cat > cpptraj.in <<EOF
parm ../tz2.nhe.parm7
ensemblesize 4
ensemble rem.crd.001 nosort
distance EndToEnd :1 :12 out distances.dat
run
runanalysis avg EndToEnd out avg.distances.dat name AvgE2E
EOF
RunCpptraj "$UNITNAME"
if [ -z "$DO_PARALLEL" ] ; then
  DoTest distances.dat.save distances.dat
  DoTest avg.distances.dat.save avg.distances.dat
else
  DoTest distances.dat.0.save distances.dat.0
  DoTest distances.dat.1.save distances.dat.1
  DoTest distances.dat.2.save distances.dat.2
  DoTest distances.dat.3.save distances.dat.3
  DoTest avg.distances.dat.0.save avg.distances.dat.0
  DoTest avg.distances.dat.1.save avg.distances.dat.1
  DoTest avg.distances.dat.2.save avg.distances.dat.2
  DoTest avg.distances.dat.3.save avg.distances.dat.3
fi

EndTest
