#!/bin/bash

. ../MasterTest.sh

CleanFiles corr.in v1.auto.dat v1v2.cross.dat v1.dplr.auto.dat \
           v1v2.dplr.cross.dat tz2.2.3.cross.dat \
           v1.o0.auto.dat v1.o1.auto.dat

TESTNAME='Timecorr test'
Requires netcdf
INPUT="corr.in"
TOP=../tz2.parm7

# Time correlation functions
cat > corr.in <<EOF
trajin ../tz2.nc
vector v1 :2  :4
vector v2 :12 :10
 
analyze timecorr vec1 v1 out v1.auto.dat ptrajformat
timecorr vec1 v1 out v1.o1.auto.dat order 1
timecorr vec1 v1 out v1.o0.auto.dat order 0
analyze timecorr vec1 v1 vec2 v2 out v1v2.cross.dat ptrajformat
analyze timecorr vec1 v1 dplr out v1.dplr.auto.dat ptrajformat
analyze timecorr vec1 v1 vec2 v2 dplr out v1v2.dplr.cross.dat ptrajformat

# Dipole correlation analysis
vector v3 @14 @15
vector v4 @38 @39
analyze timecorr vec1 v3 vec2 v4 dplr out tz2.2.3.cross.dat ptrajformat
EOF
RunCpptraj "Timecorr test."
DoTest v1.auto.dat.save v1.auto.dat
DoTest v1v2.cross.dat.save v1v2.cross.dat
DoTest v1.dplr.auto.dat.save v1.dplr.auto.dat
DoTest v1v2.dplr.cross.dat.save v1v2.dplr.cross.dat
DoTest v1.o0.auto.dat.save v1.o0.auto.dat
DoTest v1.o1.auto.dat.save v1.o1.auto.dat
DoTest tz2.2.3.cross.dat.save tz2.2.3.cross.dat

EndTest
exit 0
