#!/bin/bash

. ../MasterTest.sh

# Clean
# NOTE: CpptrajPairDist name defined in Action_Clustering.cpp
CleanFiles cluster.in cnumvtime.dat avg.summary.dat summary.dat CpptrajPairDist cpop.agr

# Test 1
CheckNetcdf
cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc
cluster C1 :2-10 clusters 3 epsilon 4.0 out cnumvtime.dat summary avg.summary.dat nofit savepairdist cpopvtime cpop.agr 
cluster crd1 :2-10 clusters 3 epsilon 4.0 summary summary.dat complete nofit loadpairdist
EOF
INPUT="-i cluster.in"
RunCpptraj "Cluster command test."
DoTest cnumvtime.dat.save cnumvtime.dat
DoTest avg.summary.dat.save avg.summary.dat 
DoTest summary.dat.save summary.dat
DoTest cpop.agr.save cpop.agr
CheckTest

EndTest

exit 0
