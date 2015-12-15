#!/bin/bash

. ../MasterTest.sh

CleanFiles cluster.in symmrmsd.*.dat rms*.dat srmsd*.dat 2drms.gnu

INPUT="-i cluster.in"

# Run Clustering
RunCluster() {
  echo "parm ../Test_SymmRmsd/AFV.parm7" > cluster.in
  for ((i=0; i < 10; i++)) ; do
    echo "trajin ../Test_SymmRmsd/AFV.rst7" >> cluster.in
  done
  for ((i=0; i < 10; i++)) ; do
    echo "trajin ../Test_SymmRmsd/AFV.rotate2.rst7" >> cluster.in
  done
  for ((i=0; i < 10; i++)) ; do
    echo "trajin ../Test_SymmRmsd/AFV.rotate3.rst7" >> cluster.in
  done
  for ((i=0; i < 10; i++)) ; do
    echo "trajin ../Test_SymmRmsd/AFV.rotate23.rst7" >> cluster.in
  done
  cat >> cluster.in <<EOF
2drms !@H= out 2drms.gnu srmsd
cluster cluster1 hieragglo epsilon 0.8 averagelinkage \
        rms !@H= \
        out rms.dat summary rms.summary.dat info rms.info.dat
#debug analysis 2
cluster cluster2 hieragglo epsilon 0.8 averagelinkage \
        srmsd !@H= \
        out srmsd.dat summary srmsd.summary.dat info srmsd.info.dat
EOF
  RunCpptraj "Clustering with symmetry-corrected RMSD metric (also 2D SRMSD)"
  DoTest rms.summary.dat.save rms.summary.dat
  DoTest rms.info.dat.save rms.info.dat
  DoTest srmsd.summary.dat.save srmsd.summary.dat
  DoTest srmsd.info.dat.save srmsd.info.dat
  DoTest 2drms.gnu.save 2drms.gnu
}

RunCluster

EndTest
exit 0
