#!/bin/bash

. ../MasterTest.sh

CleanFiles dbscan.in summary.dat info.dat rvc.dat sievesummary.dat.? \
  sieveinfo.dat.? Kdist.dat Kdist.4.dat rms2d.gnu Kmatrix.gnu Kmatrix.max.dat
INPUT="-i dbscan.in"
TESTNAME='Cluster DBSCAN tests'
Requires netcdf
# Test clustering
cat > dbscan.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
createcrd crd1
#debug analysis 1
cluster crdset crd1 C0 @CA dbscan epsilon 1.7 minpoints 5 summary summary.dat info info.dat gracecolor
#2drms crdset crd1 @CA rmsout rms2d.gnu 
rms R0 first @CA
create rvc.dat R0 C0
EOF
RunCpptraj "DBSCAN test"
DoTest rvc.dat.save rvc.dat
CheckTest

# Test 4-dist plot generation 
cat > dbscan.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
cluster C0 @CA dbscan kdist 4 #summary sievesummary.dat.1 info sieveinfo.dat.1 epsilon 1.9 minpoints 4
EOF
RunCpptraj "DBSCAN automatic parameter determination test." 
DoTest Kdist.dat.save Kdist.4.dat

# Test kdist map generation
#for ((KVAL=1 ; KVAL <= 20; KVAL++)) ; do
#  cat > dbscan.in <<EOF
#parm ../tz2.parm7
#trajin ../tz2.nc
#cluster C0 @CA dbscan kdist $KVAL
#EOF
#  RunCpptraj "DBSCAN kdist map test $KVAL."
#done
#mv Kdist.*.dat Kdist/
cat > dbscan.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
cluster C0 @CA dbscan kdist 1-20
EOF
RunCpptraj "DBSCAN kdist map test"
DoTest Kmatrix.gnu.save Kmatrix.gnu

# Test with sieving
cat > dbscan.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
cluster C0 @CA dbscan epsilon 1.7 minpoints 5 bestrep cumulative \
        summary sievesummary.dat.2 info sieveinfo.dat.2 sieve 5
EOF
RunCpptraj "DBSCAN with sieve"
DoTest sieveinfo.dat.2.save sieveinfo.dat.2
CheckTest

EndTest
exit 0
