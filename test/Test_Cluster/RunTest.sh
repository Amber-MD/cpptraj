#!/bin/bash

. ../MasterTest.sh

# Clean
# NOTE: CpptrajPairDist name defined in Action_Clustering.cpp
CleanFiles cluster.in cnumvtime.dat avg.summary.dat summary.dat CpptrajPairDist \
           cpop.agr summary2.dat Cmatrix.nccmatrix Cmatrix.cmatrix summary3.dat \
           normpop.agr normframe.agr cascii.dat.save cascii.dat pw.out \
           cinfo.dat mysil.cluster.dat mysil.frame.dat cascii?.info

TESTNAME='Hierarchical agglomerative clustering tests'
#Requires netcdf
INPUT="-i cluster.in"
# Test in-memory PW dist calc
cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.crd
cluster C1 :2-10 clusters 3 epsilon 4.0 out cnumvtime.dat info cinfo.dat sil mysil summary avg.summary.dat nofit savepairdist cpopvtime cpop.agr pairdist Cmatrix.cmatrix 
cluster crd1 :2-10 clusters 3 epsilon 4.0 summary summary.dat complete nofit loadpairdist pairdist Cmatrix.cmatrix
EOF
RunCpptraj "Cluster command test, in-memory pairwise distances."
DoTest cnumvtime.dat.save cnumvtime.dat
DoTest cinfo.dat.save cinfo.dat
DoTest mysil.cluster.dat.save mysil.cluster.dat
DoTest mysil.frame.dat.save mysil.frame.dat
DoTest avg.summary.dat.save avg.summary.dat 
DoTest summary.dat.save summary.dat
DoTest cpop.agr.save cpop.agr

# Test loading PW distances from Cmatrix file
# If netcdf support, also write out netcdf cmatrix file
if [ ! -z "$CPPTRAJ_NETCDFLIB" ] ; then
  NCWRITE='writedata Cmatrix.nccmatrix PW'
fi
cat > cluster.in <<EOF
readdata Cmatrix.cmatrix name PW
parm ../tz2.parm7
loadtraj ../tz2.crd name MyTraj
runanalysis cluster crd1 crdset MyTraj :2-10 clusters 3 epsilon 4.0 summary summary2.dat \
                    complete nofit pairdist PW \
                    cpopvtime normpop.agr normpop
$NCWRITE
EOF
RunCpptraj "Cluster command test, read pairwise distances."
DoTest summary.dat.save summary2.dat
DoTest normpop.agr.save normpop.agr

# Test loading PW distances from NetCDF cmatrix file
UNITNAME='Cluster command test, read NetCDF pairwise distances'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > cluster.in <<EOF
readdata Cmatrix.nccmatrix name PW
parm ../tz2.parm7
loadtraj ../tz2.crd name MyTraj
runanalysis cluster crd1 crdset MyTraj :2-10 clusters 3 epsilon 4.0 summary summary3.dat \
                    complete nofit pairdist PW \
                    cpopvtime normframe.agr normframe
writedata Cmatrix.nccmatrix PW
EOF
  RunCpptraj "$UNITNAME."
  DoTest summary.dat.save summary3.dat
  DoTest normframe.agr.save normframe.agr
fi

# Test writing/reading ASCII cluster pairwise file
cat > cluster.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd
cluster C1 :2-10 clusters 3 epsilon 4.0 out cascii.dat.save nofit savepairdist pairdist pw.out \
  sieve 6 random sieveseed 2 info cascii1.info
EOF
RunCpptraj "Cluster command test, write ASCII pairwise distances."
cat > cluster.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd
cluster C1 :2-10 clusters 3 epsilon 4.0 out cascii.dat nofit loadpairdist pairdist pw.out \
  info cascii2.info
EOF
RunCpptraj "Cluster command test, read ASCII pairwise distances."
DoTest cascii.dat.save cascii.dat
DoTest cascii1.info cascii2.info

EndTest

exit 0
