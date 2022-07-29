#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cluster.in C?.cnumvtime.dat C?.cinfo.dat C?.summary.dat \
           C?.cmatrix.dat

TESTNAME='Clustering with quaternion RMSD tests'

INPUT="-i cluster.in"

UNITNAME='HA clustering with regular RMSD calc.'
cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.crd
cluster C1 :2-10 rms clusters 3 epsilon 4.0 \
  out C1.cnumvtime.dat info C1.cinfo.dat summary C1.summary.dat \
  savepairdist pairdist C1.cmatrix.dat
EOF
RunCpptraj "$UNITNAME"

UNITNAME='HA clustering with quaternion RMSD calc.'
cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.crd
cluster C1 :2-10 qrmsd clusters 3 epsilon 4.0 \
  out C2.cnumvtime.dat info C2.cinfo.dat summary C2.summary.dat \
  savepairdist pairdist C2.cmatrix.dat
EOF
RunCpptraj "$UNITNAME"
DoTest C1.cnumvtime.dat C2.cnumvtime.dat
DoTest C1.cinfo.dat C2.cinfo.dat
DoTest C1.summary.dat C2.summary.dat


EndTest
