#!/bin/bash

. ../MasterTest.sh

CleanFiles image.in reimage.mdcrd image.G3_3A.rst7

TRAJ=ptraj.image.nc
INPUT="-i image.in"

TESTNAME='AutoImage test'
RequiresMaxThreads 2 "$TESTNAME"

CheckNetcdf "$TESTNAME"
if [ $? -ne 0 ] ; then
  SkipCheck "$TESTNAME"
else
  cat > image.in <<EOF
parm ../dna30.parm7
trajin split.duplex.nc
autoimage firstatom
trajout reimage.mdcrd
EOF
  RunCpptraj "AutoImage test"
  DoTest reimage.mdcrd.save reimage.mdcrd
fi

TESTNAME="AutoImage test, use anchor mask"
MaxThreads 1 "$TESTNAME"
if [ $? -ne 0 ] ; then
  SkipCheck "$TESTNAME"
else
  cat > image.in <<EOF
parm nowat.G3_3A.parm7
trajin G3_3A.rst7
autoimage anchor :96 origin
trajout image.G3_3A.rst7
EOF
  RunCpptraj "$TESTNAME"
  DoTest image.G3_3A.rst7.save image.G3_3A.rst7
fi

EndTest

exit 0 
