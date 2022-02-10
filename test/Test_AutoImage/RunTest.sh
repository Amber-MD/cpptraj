#!/bin/bash

. ../MasterTest.sh

CleanFiles image.in reimage.mdcrd image.G3_3A.rst7 image.cpptraj.Fg238.rst7

TRAJ=ptraj.image.nc
INPUT="-i image.in"

TESTNAME='AutoImage tests'
Requires maxthreads 2

UNITNAME='AutoImage test (split DNA duplex)'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > image.in <<EOF
parm ../dna30.parm7
trajin split.duplex.nc
autoimage firstatom
trajout reimage.mdcrd
EOF
  RunCpptraj "$UNITNAME"
  DoTest reimage.mdcrd.save reimage.mdcrd
fi

UNITNAME="AutoImage test, use anchor mask"
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > image.in <<EOF
parm nowat.G3_3A.parm7
trajin G3_3A.rst7
autoimage anchor :96 origin
trajout image.G3_3A.rst7
EOF
  RunCpptraj "$UNITNAME"
  DoTest image.G3_3A.rst7.save image.G3_3A.rst7
fi

UNITNAME="AutoImage test, non-sequential molecules"
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > image.in <<EOF
parm cpptraj.Fg238.parm7
trajin cpptraj.Fg238.rst7
autoimage
trajout image.cpptraj.Fg238.rst7
EOF
  RunCpptraj "$UNITNAME"
  DoTest image.cpptraj.Fg238.rst7.save image.cpptraj.Fg238.rst7
fi

EndTest

exit 0 
