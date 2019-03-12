#!/bin/bash

. ../MasterTest.sh

in=density.in
out1=number_density.dat
out2=mass_density.dat
out3=charge_density.dat
out4=electron_density.dat

CleanFiles $in $out1 $out2 $out3 $out4 total.dat tz2.wato.agr

INPUT="-i $in"
TESTNAME='Density tests'
Requires maxthreads 10

# Density along an axis
UNITNAME='Density along axis tests'
CheckFor maxthreads 1 
if [ $? -eq 0 ] ; then
  del='delta 0.25'
  masks='":PC@P31" ":PC@N31" ":PC@C2" ":PC | :OL | :OL2"'

  cat > $in <<EOF
# crd/top courtesy of Callum Dickson, Imperial College London
parm ../DOPC.parm7
trajin ../DOPC.rst7

center ":PC | :OL | :OL2" origin

density out $out1 number $del $masks
density out $out2 mass $del $masks
density out $out3 charge $del $masks
density out $out4 electron $del $masks
EOF
  RunCpptraj "$UNITNAME"
  DoTest ${out1}.save $out1
  DoTest ${out2}.save $out2
  DoTest ${out3}.save $out3
  DoTest ${out4}.save $out4
fi

# Density along an axis, multiple frames
UNITNAME='Multi-frame density along axis test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
   cat > $in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
autoimage origin
density :WAT@O out tz2.wato.agr xydy
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2.wato.agr.save tz2.wato.agr
fi

# Total system density test
UNITNAME='Total system density test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > $in <<EOF
parm ../tz2.truncoct.parm7 [OCT]
trajin ../tz2.truncoct.nc parm [OCT]
density D1
go
clear trajin

parm ../tz2.ortho.parm7 [ORTHO]
trajin ../tz2.ortho.nc parm [ORTHO]
density D2
go

writedata total.dat D1 D2
EOF
  RunCpptraj "$UNITNAME"
  DoTest total.dat.save total.dat
fi

EndTest

exit 0
