#!/bin/bash

. ../MasterTest.sh

CleanFiles cluster.in symmrmsd.*.dat rms*.dat srmsd*.dat 2drms.gnu

INPUT="-i cluster.in"
TESTNAME='Clustering with symmetry-corrected RMSD metric'
Requires netcdf
# Run Clustering
RunCluster() {
  cat > cluster.in <<EOF
parm ../AFV.parm7
trajin ../AFV.nc
2drms !@H= out 2drms.gnu srmsd
cluster cluster1 hieragglo epsilon 0.8 averagelinkage \
        rms !@H= \
        out rms.dat summary rms.summary.dat info rms.info.dat
#debug analysis 2
cluster cluster2 hieragglo epsilon 0.8 averagelinkage \
        srmsd !@H= \
        out srmsd.dat summary srmsd.summary.dat info srmsd.info.dat
EOF
  RunCpptraj "$TESTNAME (also 2D SRMSD)"
  DoTest rms.summary.dat.save rms.summary.dat
  DoTest rms.info.dat.save rms.info.dat
  DoTest srmsd.summary.dat.save srmsd.summary.dat
  DoTest srmsd.info.dat.save srmsd.info.dat
  DoTest 2drms.gnu.save 2drms.gnu
}

RunCluster

EndTest
exit 0
