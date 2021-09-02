#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in flatten.dat

INPUT='-i cpptraj.in'

cat > cpptraj.in <<EOF
readdata ../Test_Matrix/mtest.dat.13.save name Mat read2d
flatten Mat name FlatSum mode sum
flatten Mat name FlatAvg mode avg
writedata flatten.dat FlatSum FlatAvg
EOF
RunCpptraj "Flatten test."
DoTest flatten.dat.save flatten.dat

EndTest
