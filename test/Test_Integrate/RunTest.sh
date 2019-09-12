#!/bin/bash

. ../MasterTest.sh

CleanFiles curve.in curve.dat int.dat

INPUT="-i curve.in"

cat > curve.in <<EOF
readdata ../Test_Corr/corr.dat.save name Corr
integrate Corr out curve.dat name Int_Corr intout int.dat
datafile int.dat prec 16.8
EOF
RunCpptraj "Integration test."
DoTest curve.dat.save curve.dat

EndTest
exit 0
