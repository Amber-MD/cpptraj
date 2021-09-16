#!/bin/bash

. ../MasterTest.sh

TESTNAME='Test clustering on coordinates and data'

CleanFiles cluster.in c1.info.dat c1.summary.dat c1.cnvt.dat data.dat metricstats.dat

INPUT='-i cluster.in'

cat > cluster.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.crd name TZ2

crdaction TZ2 distance DATA :1 :13 out data.dat

runanalysis cluster c1 \
  hieragglo clusters 5 \
  data TZ2,DATA \
  rms :1-12&!@H= \
  out c1.cnvt.dat \
  info c1.info.dat \
  summary c1.summary.dat metricstats metricstats.dat
EOF
RunCpptraj "$TESTNAME"
DoTest c1.info.dat.save c1.info.dat
DoTest c1.summary.dat.save c1.summary.dat
DoTest c1.cnvt.dat.save c1.cnvt.dat
DoTest metricstats.dat.save metricstats.dat

EndTest
