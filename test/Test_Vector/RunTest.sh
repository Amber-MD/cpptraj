#!/bin/bash

. ../MasterTest.sh

CleanFiles vector.in vtest.dat.? v8.mol2 corr.v0.v8.dat

CheckNetcdf

INPUT="-i vector.in"
# Test Vector mask, principle xyz, dipole, box 
cat > vector.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
vector v0 principal x @CA out vtest.dat.0 ptrajoutput
vector v1 principal y @CA out vtest.dat.1 ptrajoutput
vector v2 principal z @CA out vtest.dat.2 ptrajoutput
vector v3 @91 @92 out vtest.dat.3 ptrajoutput
vector v4 @91 dipole out vtest.dat.4 ptrajoutput
vector v5 box out vtest.dat.5 ptrajoutput
vector v6 center out vtest.dat.6 :1
vector v7 corrplane out vtest.dat.7 :2@CD2,CE?,CZ?,CH2
vector v8 minimage out vtest.dat.8 :4 :11
corr v0 v8 out corr.v0.v8.dat
run
writedata v8.mol2 vectraj v8 trajfmt mol2
EOF
RunCpptraj "Vector Tests"
DoTest vtest.dat.0.save vtest.dat.0
DoTest vtest.dat.1.save vtest.dat.1
DoTest vtest.dat.2.save vtest.dat.2
DoTest vtest.dat.3.save vtest.dat.3
DoTest vtest.dat.4.save vtest.dat.4
DoTest vtest.dat.5.save vtest.dat.5
DoTest vtest.dat.6.save vtest.dat.6
DoTest vtest.dat.7.save vtest.dat.7
DoTest v8.mol2.save v8.mol2
CheckTest

EndTest
  
exit 0
