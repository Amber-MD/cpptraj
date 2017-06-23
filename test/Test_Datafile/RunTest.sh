#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles prec.in prec.dat a1.dat a1.agr xprec.dat

CheckNetcdf
TOP="../tz2.truncoct.parm7"

# Test 1
cat > prec.in <<EOF
noprogress
trajin ../tz2.truncoct.nc
rms R0 first :2-11 out prec.dat
precision prec.dat * 8 3
EOF
INPUT="prec.in"
RunCpptraj "Data file output precision test."
DoTest prec.dat.save prec.dat
CheckTest

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

EndTest

exit 0
