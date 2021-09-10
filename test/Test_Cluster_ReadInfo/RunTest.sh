#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cluster.in C?.info.dat PW

TESTNAME='Cluster read info tests'
#Requires netcdf
INPUT="-i cluster.in"

cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
loadcrd ../tz2.crd name MyCrd

# Cluster to 3 clusters
cluster crdset MyCrd C1 :2-10 hieragglo clusters 3 info C1.info.dat savepairdist pairdist PW
# Cluster to 8 clusters, use existing PW
cluster crdset MyCrd C2 :2-10 hieragglo clusters 8 info C2.info.dat pairdist PW
EOF
RunCpptraj "Cluster restart test, part 1"

cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
loadcrd ../tz2.crd name MyCrd

# Cluster to 3 clusters, start from existing info with 8 clusters
cluster crdset MyCrd C3 :2-10 hieragglo clusters 3 \
  readinfo infofile C2.info.dat \
  info C3.info.dat \
  loadpairdist pairdist PW
EOF
RunCpptraj "Cluster restart test, part 2"
DoTest C1.info.dat C3.info.dat

EndTest

exit 0
