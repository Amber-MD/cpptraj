#!/bin/bash

. ../MasterTest.sh

CleanFiles vectors.dat dotproduct.dat corr.in v1init_dot_v1.dat Magnitude.dat

INPUT="corr.in"
TOP=../tz2.parm7

TESTNAME='Vector math tests'

# Dot/cross products
UNITNAME='Vector dot/cross product tests'
Requires netcdf
if [ $? -eq 0 ] ; then
  cat > corr.in <<EOF
readdata v1initial.dat vector
trajin ../tz2.nc
vector v1 :2  :4  #out vectors.dat
vector v2 :12 :10 #out vectors.dat
vectormath vec1 v1 vec2 v2 dotproduct out dotproduct.dat name V1*V2
vectormath vec1 v1 vec2 v2 dotproduct norm out dotproduct.dat name |V1|*|V2|
vectormath vec1 v1 vec2 v2 dotangle out dotproduct.dat name acos(|V1|*|V2|)
vectormath vec1 v1 vec2 v2 crossproduct out dotproduct.dat name |V1|x|V2|
vectormath vec1 v1initial.dat vec2 v1 dotangle out v1init_dot_v1.dat name V1o_V1_ang
EOF
  RunCpptraj "$UNITNAME"
  DoTest dotproduct.dat.save dotproduct.dat
  DoTest v1init_dot_v1.dat.save v1init_dot_v1.dat
fi

UNITNAME='Vector magnitude tests'
cat > corr.in <<EOF
readdata ../Test_Vector/avgcoord.out.save vector name V1 magnitude
runanalysis vectormath vec1 V1 magnitude out Magnitude.dat name V1Mag
EOF
RunCpptraj "$UNITNAME"
DoTest Magnitude.dat.save Magnitude.dat

EndTest
exit 0 
