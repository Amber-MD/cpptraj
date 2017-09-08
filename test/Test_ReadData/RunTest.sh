#!/bin/bash

. ../MasterTest.sh

CleanFiles vector.in v6and7.dat

INPUT="-i vector.in"
# Test read/append of vector dataset 
cat > vector.in <<EOF
readdata ../Test_Vector/vtest.dat.6.save vector name v6and7
readdata ../Test_Vector/vtest.dat.7.save vector name v6and7
writedata v6and7.dat v6and7
EOF
RunCpptraj "Read vector data test"
DoTest v6and7.dat.save v6and7.dat

EndTest
  
exit 0
