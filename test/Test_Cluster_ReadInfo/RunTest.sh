#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cluster.in C?.info.dat PW C2.dat C??.info.dat

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
cluster crdset MyCrd C2 :2-10 hieragglo clusters 8 info C2.info.dat pairdist PW out C2.dat
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

# Same but use cnumvime set
readdata C2.dat name MyClusters
cluster crdset MyCrd C10 :2-10 hieragglo clusters 3 \
  readinfo cnvtset MyClusters \
  info C10.info.dat \
  pairdist PW
EOF
RunCpptraj "Hierarchical agglomerative cluster restart test, part 2"
DoTest C1.info.dat C3.info.dat
DoTest C1.info.dat C10.info.dat

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
#debug 10
# Cluster to 3 clusters, start from existing info with 8 clusters
cluster crdset MyCrd C6 :2-10 kmeans clusters 3 \
  readinfo infofile C5.info.dat \
  info C6.info.dat \
  loadpairdist pairdist PW
EOF
RunCpptraj "K-means cluster restart test, part 2"
DoTest C4.info.dat C6.info.dat

# K-means, try to cluster to more clusters from fewer clusters.
cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
loadcrd ../tz2.crd name MyCrd

cluster crdset MyCrd C7 :2-10 kmeans clusters 8 \
  readinfo infofile C4.info.dat \
  info C7.info.dat \
  loadpairdist pairdist PW
EOF
RunCpptraj "K-means cluster restart, 3 to 8 clusters."
DoTest C7.info.dat.save C7.info.dat

# DBscan
cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
loadcrd ../tz2.crd name MyCrd

cluster crdset MyCrd C8 @CA dbscan epsilon 1.0 minpoints 5 info C8.info.dat \
  pairdist PW.CA savepairdist
EOF
RunCpptraj "DBscan cluster restart, part 1"

cat > cluster.in <<EOF
noprogress
parm ../tz2.parm7
loadcrd ../tz2.crd name MyCrd

cluster crdset MyCrd C9 @CA dbscan epsilon 1.7 minpoints 5 info C9.info.dat \
  readinfo infofile C8.info.dat \
  pairdist PW.CA loadpairdist
EOF
RunCpptraj "DBscan cluster restart, part 2"
DoTest C9.info.dat.save C9.info.dat

EndTest

exit 0
