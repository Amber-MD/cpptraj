#!/bin/bash

. ../MasterTest.sh

CleanFiles lowest.in lowest.dat All.agr

INPUT="-i lowest.in"

Test1() {
cat > lowest.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
rms myrmsd @CA
lowestcurve points 10 myrmsd out lowest.dat
EOF
RunCpptraj "LowestCurve test"
}

cat > lowest.in <<EOF
readdata esurf_vs_rmsd.dat.txt index 1 name MyData
list dataset
runanalysis lowestcurve MyData points 10 step 0.2 name Lowest
writedata All.agr MyData Lowest
EOF
RunCpptraj "LowestCurve test"
DoTest All.agr.save All.agr

EndTest
exit 0
