#!/bin/bash

. ../MasterTest.sh

# Clean
# NOTE: CpptrajPairDist name defined in Action_Clustering.cpp
CleanFiles cluster.in cnumvtime.dat avg.summary.dat summary.dat CpptrajPairDist \
           cpop.agr summary2.dat Cmatrix.nccmatrix Cmatrix.cmatrix summary3.dat \
           normpop.agr normframe.agr

TESTNAME='Hierarchical agglomerative clustering tests'
Requires netcdf
INPUT="-i cluster.in"
# Test in-memory PW dist calc
cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc
cluster C1 :2-10 clusters 3 epsilon 4.0 out cnumvtime.dat summary avg.summary.dat nofit savepairdist cpopvtime cpop.agr pairdist Cmatrix.cmatrix 
cluster crd1 :2-10 clusters 3 epsilon 4.0 summary summary.dat complete nofit loadpairdist pairdist Cmatrix.cmatrix
EOF
RunCpptraj "Cluster command test, in-memory pairwise distances."
DoTest cnumvtime.dat.save cnumvtime.dat
DoTest avg.summary.dat.save avg.summary.dat 
DoTest summary.dat.save summary.dat
DoTest cpop.agr.save cpop.agr
# Test loading PW distances from Cmatrix file
cat > cluster.in <<EOF
readdata Cmatrix.cmatrix name PW
parm ../tz2.parm7
loadtraj ../tz2.nc name MyTraj
runanalysis cluster crd1 crdset MyTraj :2-10 clusters 3 epsilon 4.0 summary summary2.dat \
                    complete nofit pairdist PW \
                    cpopvtime normpop.agr normpop
writedata Cmatrix.nccmatrix PW
EOF
RunCpptraj "Cluster command test, read pairwise distances."
DoTest summary.dat.save summary2.dat
DoTest normpop.agr.save normpop.agr
# Test loading PW distances from NetCDF cmatrix file
cat > cluster.in <<EOF
readdata Cmatrix.nccmatrix name PW
parm ../tz2.parm7
loadtraj ../tz2.nc name MyTraj
runanalysis cluster crd1 crdset MyTraj :2-10 clusters 3 epsilon 4.0 summary summary3.dat \
                    complete nofit pairdist PW \
                    cpopvtime normframe.agr normframe
writedata Cmatrix.nccmatrix PW
EOF
RunCpptraj "Cluster command test, read NetCDF pairwise distances."
DoTest summary.dat.save summary3.dat
DoTest normframe.agr.save normframe.agr

EndTest

exit 0
