#!/bin/bash

. ../MasterTest.sh

CleanFiles spline.in spline.dat

INPUT="-i spline.in"

cat > spline.in <<EOF
readdata ../Test_Temperature/T2.dat.save name T2
spline T2 out spline.dat meshsize 100 meshmin 0.0 meshmax 11.0 
EOF
RunCpptraj "Spline test"
DoTest spline.dat.save spline.dat
EndTest
exit 0
