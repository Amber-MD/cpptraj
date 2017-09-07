#!/bin/bash

. ../MasterTest.sh

CleanFiles check.in report.dat nprob.dat tz2.dat skip.dat around.dat

INPUT="-i check.in"

RequiresMaxThreads 10 "Structure check tests"

MaxThreads 1 "Structure Check"
if [ $? -ne 0 ] ; then
  SkipCheck "Structure Check"
else
  # Test 1
  cat > check.in <<EOF
parm ../tz2.parm7
trajin tz2.stretched.pdb
check reportfile report.dat offset 0.7 out nprob.dat Tz2Check
EOF
  RunCpptraj "Structure Check"
  DoTest report.dat.save report.dat
  DoTest nprob.dat.save nprob.dat
  # Around test with skip
  TESTNAME='Structure Check with Around and Skip'
  CheckNetcdf "$TESTNAME"
  if [ $? -ne 0 ] ; then
    SkipCheck "$TESTNAME"
  else
    cat > check.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
scale :1 x 2.0 y 1.2 z 1.2
check nobondcheck :WAT around :1 out skip.dat skipbadframes silent
distance d1 out d1.dat :1 :12
EOF
    RunCpptraj "Structure Check with Around and Skip"
    DoTest skip.dat.save skip.dat
    DoTest d1.dat.save d1.dat
  fi
fi

TESTNAME='Structure check with PBC/around tests'
CheckNetcdf "$TESTNAME"
if [ $? -ne 0 ] ; then
  SkipCheck "$TESTNAME"
else
  cat > check.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
rms first :1-13
scale :1 x 2.0 y 1.2 z 1.2
check reportfile tz2.dat offset 1.0
EOF
  RunCpptraj "Structure Check with PBC"
  DoTest tz2.dat.save tz2.dat

  # Around test
  cat > check.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
rms first :1-13
scale :1 x 2.0 y 1.2 z 1.2
check reportfile around.dat offset 1.0 :WAT around :1
EOF
  RunCpptraj "Structure Check with Around"
  DoTest around.dat.save around.dat
fi
EndTest

exit 0
