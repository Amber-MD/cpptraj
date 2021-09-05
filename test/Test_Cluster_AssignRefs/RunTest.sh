#!/bin/bash

. ../MasterTest.sh

CleanFiles cluster.in summary.dat ksummary.dat

INPUT='-i cluster.in'

TESTNAME='Clustering, assign reference names'

cat > cluster.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc

cluster Tz2 rms :1-12@N,CA,C,O hieragglo clusters 5 summary summary.dat pairdist PW
cluster Kmeans rms :1-12@N,CA,C,O kmeans clusters 5 summary ksummary.dat pairdist PW
EOF
RunCpptraj "Test1"


EndTest
