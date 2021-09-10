#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cluster.in C?.info.dat PW

TESTNAME='Cluster read info tests'
#Requires netcdf
INPUT="-i cluster.in"

# Hier. Agglo.
cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
loadcrd ../tz2.crd name MyCrd

# Cluster to 3 clusters
cluster crdset MyCrd C1 :2-10 hieragglo clusters 3 info C1.info.dat savepairdist pairdist PW
# Cluster to 8 clusters, use existing PW
cluster crdset MyCrd C2 :2-10 hieragglo clusters 8 info C2.info.dat pairdist PW
EOF
RunCpptraj "Hierarchical aggolmerative cluster restart test, part 1"

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
RunCpptraj "Hierarchical agglomerative cluster restart test, part 2"
DoTest C1.info.dat C3.info.dat

# K-means
cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
loadcrd ../tz2.crd name MyCrd

# Cluster to 3 clusters
cluster crdset MyCrd C4 :2-10 kmeans clusters 3 info C4.info.dat savepairdist pairdist PW
# Cluster to 8 clusters, use existing PW
cluster crdset MyCrd C5 :2-10 kmeans clusters 8 info C5.info.dat pairdist PW
EOF
RunCpptraj "K-means cluster restart test, part 1"

cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
loadcrd ../tz2.crd name MyCrd
debug 10
# Cluster to 3 clusters, start from existing info with 8 clusters
cluster crdset MyCrd C6 :2-10 kmeans clusters 3 \
  readinfo infofile C5.info.dat \
  info C6.info.dat \
  loadpairdist pairdist PW
EOF
RunCpptraj "K-means cluster restart test, part 2"
DoTest C4.info.dat C6.info.dat

EndTest

exit 0
