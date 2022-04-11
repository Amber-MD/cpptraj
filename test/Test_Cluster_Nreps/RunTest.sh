#!/bin/bash

. ../MasterTest.sh

CleanFiles cluster.in 2drms.gnu clusters.*.dat summary.*.dat singlerep.nc \
           Rep.c*.nc PD cumulative centroid cumulative_nosieve

TESTNAME='Cluster representative frames tests'
INPUT='-i cluster.in'

UNITNAME='Test saving 5 best reps'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > cluster.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
#2drms !@H= out 2drms.gnu
cluster kmeans clusters 2 info clusters.2.dat summary summary.2.dat \
        savepairdist pairdist PD savenreps 5 
EOF
  RunCpptraj "$TESTNAME"

  DoTest clusters.2.dat.save clusters.2.dat
  DoTest summary.2.dat.save summary.2.dat
fi

UNITNAME='Test cumulative best rep'
cat > cluster.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd
cluster hieragglo clusters 5 rms @CA bestrep cumulative singlerepout cumulative
EOF
RunCpptraj "$UNITNAME"
DoTest cumulative.save cumulative

UNITNAME='Test centroid best rep'
cat > cluster.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd
cluster hieragglo clusters 5 rms @CA bestrep centroid singlerepout centroid
EOF
RunCpptraj "$UNITNAME"
DoTest centroid.save centroid

UNITNAME='Test cumulative, no sieve best rep'
cat > cluster.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd
cluster hieragglo clusters 5 rms @CA sieve 5 bestrep cumulative_nosieve singlerepout cumulative_nosieve
EOF
RunCpptraj "$UNITNAME"
DoTest cumulative_nosieve.save cumulative_nosieve

EndTest
exit 0
