#!/bin/bash

. ../MasterTest.sh

CleanFiles box.in addbox.rst7 addbox.rst7.? addbox.rst7.10 \
                  modX.rst7   modX.rst7.?   modX.rst7.10 \
                  frame1.rst7 tz2.box.rst7 tz2.vdw.rst7 \
                  *.ucell.dat *.frac.dat *.shape.dat

TESTNAME='Box tests'
Requires netcdf maxthreads 10

INPUT="-i box.in"

UNITNAME='Box Test (Add box info)'
cat > box.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10 
strip !(:1)
box x 42.428  y 42.428  z 42.428 alpha 109.471 beta 109.471 gamma 109.471
trajout addbox.rst7 time0 0
go
EOF
RunCpptraj "$UNITNAME"
DoTest addbox.rst7.1.save addbox.rst7.1
DoTest addbox.rst7.10.save addbox.rst7.10

UNITNAME='Box test (remove box info)'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > box.in <<EOF
parm ../tz2.parm7
parmstrip !(:1)
trajin addbox.rst7.1.save
box nobox
trajout frame1.rst7
EOF
  RunCpptraj "$UNITNAME"
  DoTest frame1.rst7.save frame1.rst7
fi

UNITNAME='Box test (Modify box length)'
cat > box.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 10
strip !(:1)
box x 45.0
trajout modX.rst7 time0 0
go
EOF
RunCpptraj "$UNITNAME"
DoTest modX.rst7.1.save modX.rst7.1
DoTest modX.rst7.10.save modX.rst7.10

UNITNAME='Box test (auto orthogonal box, no radii)'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > box.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 1
box auto offset 3.0 radii none
trajout tz2.box.rst7 time0 0
run
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2.box.rst7.save tz2.box.rst7
fi

UNITNAME='Box test (auto orthogonal box, VDW radii)'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > box.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 1
box auto radii vdw
trajout tz2.vdw.rst7 time0 0
run
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2.vdw.rst7.save tz2.vdw.rst7
fi

UNITNAME='Box test (get box info)'
cat > box.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
box getbox ucell name UC out ortho.ucell.dat
box getbox frac  name FC out ortho.frac.dat
box getbox shape name SP out ortho.shape.dat
run
clear trajin
clear parm
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
box getbox ucell name UC1 out truncoct.ucell.dat
box getbox frac  name FC1 out truncoct.frac.dat
box getbox shape name SP1 out truncoct.shape.dat
EOF
RunCpptraj "$UNITNAME"
DoTest ortho.ucell.dat.save ortho.ucell.dat
DoTest ortho.frac.dat.save ortho.frac.dat
DoTest ortho.shape.dat.save ortho.shape.dat 
DoTest truncoct.ucell.dat.save truncoct.ucell.dat
DoTest truncoct.frac.dat.save truncoct.frac.dat
DoTest truncoct.shape.dat.save truncoct.shape.dat

EndTest

exit 0
