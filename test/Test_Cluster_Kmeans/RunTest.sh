#!/bin/bash

. ../MasterTest.sh

CleanFiles kmeans.in ptraj.txt ptraj.c? summary.dat info.dat cpptraj.crd.c? random.dat rinfo.dat random.crd.c?

TOP=../tz2.parm7
INPUT="kmeans.in"

RunPtraj() {
  SAVECPPTRAJ=$CPPTRAJ
  #CPPTRAJ=`which ptraj`
  CPPTRAJ=/home/droe/Amber/GIT/amber/AmberTools/src/ptraj/ptraj
  cat > kmeans.in <<EOF
trajin ../tz2.nc
cluster out ptraj means clusters 5 rms @CA verbose 1
EOF
  RunCpptraj "Ptraj Kmeans"
  DoTest ptraj.txt.save ptraj.txt
  DoTest ptraj.c0.save ptraj.c0
  mv test.out ptraj.out
}

KmeansTest() {
  cat > kmeans.in <<EOF
trajin ../tz2.nc
cluster means clusters 5 rms @CA summary summary.dat info info.dat clusterout cpptraj.crd
EOF
  RunCpptraj "Cpptraj Kmeans"
  DoTest summary.dat.save summary.dat
  DoTest info.dat.save info.dat
  DoTest cpptraj.crd.c0.save cpptraj.crd.c0
  cat > kmeans.in <<EOF
trajin ../tz2.nc
cluster means randompoint kseed 1 clusters 5 rms @CA summary random.dat info rinfo.dat clusterout random.crd
EOF
  RunCpptraj "Cpptraj Kmeans random"
  DoTest random.dat.save random.dat
  DoTest info.dat.save rinfo.dat
  DoTest cpptraj.crd.c0.save random.crd.c0
}

# ------------------------------------------------------------------------------

#RunPtraj

KmeansTest

EndTest

exit 0
