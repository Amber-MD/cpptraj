#!/bin/bash

. ../MasterTest.sh

CleanFiles vector.in vtest.dat.? vtest.dat.?? v8.mol2 corr.v0.v8.dat \
           avgcoord.out res5.out

INPUT="-i vector.in"
# Test Vector mask, principle xyz, dipole, box 
UNITNAME='Basic vector tests'
CheckFor netcdf maxthreads 10
if [ $? -eq 0 ] ; then
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
  RunCpptraj "$UNITNAME"
  DoTest vtest.dat.0.save vtest.dat.0
  DoTest vtest.dat.1.save vtest.dat.1
  DoTest vtest.dat.2.save vtest.dat.2
  DoTest vtest.dat.3.save vtest.dat.3
  DoTest vtest.dat.4.save vtest.dat.4
  DoTest vtest.dat.5.save vtest.dat.5
  DoTest vtest.dat.6.save vtest.dat.6
  DoTest vtest.dat.7.save vtest.dat.7
  DoTest v8.mol2.save v8.mol2
fi

# Test vector center with magnitude
UNITNAME='Test vector center with magnitude'
CheckFor netcdf maxthreads 10
if [ $? -eq 0 ] ; then
  cat > vector.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
vector VC_A1 center @1 out avgcoord.out magnitude
vector VC_R5 center :5 out res5.out magnitude
EOF
  RunCpptraj "$UNITNAME"
  DoTest avgcoord.out.save avgcoord.out
  DoTest res5.out.save res5.out
fi

# Test momentum vector.
UNITNAME="Momentum/velocity vector test."
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > vector.in <<EOF
parm ../tz2.parm7
trajin ../Test_SetVelocity/tz2.vel.rst7.save
vector v9 momentum out vtest.dat.9
vector v10 velocity @1-3 out vtest.dat.10
EOF
  RunCpptraj "$UNITNAME"
  DoTest vtest.dat.9.save vtest.dat.9
  DoTest vtest.dat.10.save vtest.dat.10
fi

EndTest
  
exit 0
