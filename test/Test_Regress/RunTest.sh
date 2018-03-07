#!/bin/bash

. ../MasterTest.sh

CleanFiles regress.in fit.dat statsout.dat Y.dat regress.dat

TESTNAME='Linear regression test.'

Requires maxthreads 1 

INPUT="-i regress.in"

cat > regress.in <<EOF
readdata ../Test_LowestCurve/esurf_vs_rmsd.dat.txt index 1 name XY
regress XY out regress.dat name FitXY nx 100 statsout statsout.dat
runanalysis
list dataset
createset "Y = FitXY[slope] * X + FitXY[intercept]" xstep .2 nx 100
writedata Y.dat Y
EOF
RunCpptraj "$TESTNAME"
DoTest statsout.dat.save statsout.dat
DoTest regress.dat.save regress.dat
DoTest Y.dat.save Y.dat

EndTest
exit 0
