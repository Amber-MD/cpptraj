#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles prec.in prec.dat a1.dat a1.agr xprec.dat byname.dat

TESTNAME='Data file tests'

TOP="../tz2.truncoct.parm7"
INPUT="prec.in"

# Test 1
UNITNAME='Data file output precision test'
CheckFor netcdf maxthreads 10
if [ $? -eq 0 ] ; then
  cat > prec.in <<EOF
noprogress
trajin ../tz2.truncoct.nc
rms R0 first :2-11 out prec.dat
precision prec.dat * 8 3
EOF
  RunCpptraj "$UNITNAME"
  DoTest prec.dat.save prec.dat
fi

# ReadData test
cat > prec.in <<EOF
readdata ../Test_General/a1.dat.save 
create a1.agr a1.dat.save
writedata
quit
EOF
RunCpptraj "Standard -> Grace Data"
cat > prec.in <<EOF
readdata a1.agr
write a1.dat a1.agr
quit
EOF
RunCpptraj "Grace -> Standard Data"
DoTest ../Test_General/a1.dat.save a1.dat

# xprec/xfmt
cat > prec.in <<EOF
readdata ../Test_General/a1.dat.save name A1
writedata xprec.dat A1 xprec 16.7 xfmt scientific
EOF
RunCpptraj "X column format/precision test."
DoTest xprec.dat.save xprec.dat

cat > prec.in <<EOF
readdata ../Test_RemdTraj/d1.offset.dat.save name d1
readdata ../Test_Diffusion/diff.2.dat.save index 1 name Diff
list dataset
writedata byname.dat d1 Diff groupbyname
EOF
RunCpptraj "Data file group by name test"
DoTest byname.dat.save byname.dat

EndTest

exit 0
