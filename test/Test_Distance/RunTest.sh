#!/bin/bash

. ../MasterTest.sh

CleanFiles dist.in dist.dat ortho.dat truncoct.dat \
           ortho.553.dat truncoct.385.dat

TESTNAME='Distance tests'
Requires netcdf
INPUT='-i dist.in'

cat > dist.in <<EOF
parm ../tz2.parm7
reference ../tz2.pdb
trajin ../tz2.nc

distance EndToEnd :1 :13 out dist.dat
distance ToRef @1 @1 out dist.dat reference
distance Point :1 point 0.0 0.0 0.0 out dist.dat
run
EOF
RunCpptraj "$TESTNAME"
DoTest dist.dat.save dist.dat

UNITNAME='Orthogonal imaged distance'
CheckFor maxthreads 10
if [ $? -eq 0 ] ; then
  cat > dist.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
distance EndToEnd :1 :13 out ortho.dat
distance To553         ^1 :553 out ortho.553.dat
distance To553_noImage ^1 :553 out ortho.553.dat noimage
EOF
  RunCpptraj "$UNITNAME"
  DoTest ortho.dat.save ortho.dat
  DoTest ortho.553.dat.save ortho.553.dat
fi

UNITNAME='Non-orthogonal imaged distance'
CheckFor maxthreads 10
if [ $? -eq 0 ] ; then
  cat > dist.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
distance EndToEnd :1 :13 out truncoct.dat
set RES = 385
distance To\$RES         ^1 :\$RES out truncoct.\$RES.dat
distance To\$RES_NoImage ^1 :\$RES out truncoct.\$RES.dat noimage
EOF
  RunCpptraj "$UNITNAME"
  DoTest truncoct.dat.save truncoct.dat
  DoTest truncoct.385.dat.save truncoct.385.dat
fi

EndTest
exit 0
