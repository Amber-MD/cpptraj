#!/bin/bash

. ../MasterTest.sh

CleanFiles cluster.in summary.dat ksummary.dat PW

INPUT='-i cluster.in'

TESTNAME='Clustering, assign reference names'

Requires netcdf

cat > cluster.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc

#cluster Tz2 rms :1-12@N,CA,C,O hieragglo clusters 5 summary summary.dat repout rep repfmt mol2
distance D1 :1 :12

for i=0;i<5;i++
  reference rep.c\$i.mol2 [C\$i]
done
cluster Kmeans data D1 kmeans clusters 5 summary ksummary.dat assignrefs refcut 3.0 refmask :1-12&!@H=
EOF
RunCpptraj "$TESTNAME"
DoTest ksummary.dat.save ksummary.dat

EndTest
