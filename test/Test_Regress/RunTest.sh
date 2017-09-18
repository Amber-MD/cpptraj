#!/bin/bash

. ../MasterTest.sh

CleanFiles regress.in fit.dat statsout.dat

INPUT="-i regress.in"

cat > regress.in <<EOF
readdata ../Test_LowestCurve/esurf_vs_rmsd.dat.txt index 1 name XY
regress XY out fit.dat name FitXY statsout statsout.dat
EOF
RunCpptraj "Linear regression test."
DoTest fit.dat.save fit.dat
DoTest statsout.dat.save statsout.dat

EndTest
exit 0
