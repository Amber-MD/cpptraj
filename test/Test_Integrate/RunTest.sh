#!/bin/bash

. ../MasterTest.sh

CleanFiles curve.in curve.dat

INPUT="-i curve.in"

cat > curve.in <<EOF
readdata ../Test_Corr/corr.dat.save name Corr
runanalysis integrate Corr out curve.dat name Int_Corr
EOF
RunCpptraj "Integration test."
DoTest curve.dat.save curve.dat

EndTest
exit 0
