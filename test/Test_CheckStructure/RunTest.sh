#!/bin/bash

. ../MasterTest.sh

CleanFiles check.in report.dat nprob.dat tz2.dat skip.dat around.dat d1.dat \
           partial.report.dat tz2.cut0.9.dat entangled.dat report.opcbox.dat

INPUT="-i check.in"

TESTNAME='Structure check tests'
Requires maxthreads 10

UNITNAME='Basic structure check'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  # Test 1
  cat > check.in <<EOF
parm ../tz2.parm7
trajin tz2.stretched.pdb
check reportfile report.dat offset 0.7 out nprob.dat Tz2Check
check reportfile partial.report.dat :4-6,12-13 offset 0.7
EOF
  RunCpptraj "Structure Check"
  DoTest report.dat.save report.dat
  DoTest report.dat.save partial.report.dat
  DoTest nprob.dat.save nprob.dat
  # Around test with skip
  UNITNAME='Structure Check with Around and Skip'
  CheckFor netcdf
  if [ $? -eq 0 ] ; then
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

UNITNAME='Structure check with PBC/around tests'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > check.in <<EOF
parm ../tz2.truncoct.parm7
#trajin ../tz2.truncoct.nc
trajin tz2.truncoct.baducell.nc
#rms first :1-13
#scale :1 x 2.0 y 1.2 z 1.2
check reportfile tz2.dat offset 1.0
check reportfile tz2.cut0.9.dat offset 1.0 cut 0.9
EOF
  RunCpptraj "Structure Check with PBC"
  DoTest tz2.dat.save tz2.dat
  DoTest tz2.cut0.9.dat.save tz2.cut0.9.dat
  # Around test
  cat > check.in <<EOF
parm ../tz2.truncoct.parm7
#trajin ../tz2.truncoct.nc
trajin tz2.truncoct.baducell.nc
#rms first :1-13
#scale :1 x 2.0 y 1.2 z 1.2
check reportfile around.dat offset 1.0 :WAT around :1
EOF
  RunCpptraj "Structure Check with Around"
  DoTest around.dat.save around.dat
fi

UNITNAME='Check for entangled rings'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > check.in <<EOF
parm EntangledRings.mol2
trajin EntangledRings.mol2
check reportfile entangled.dat
run
EOF
  RunCpptraj "$UNITNAME"
  DoTest entangled.dat.save entangled.dat
fi

UNITNAME='Check with extra points'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > check.in <<EOF
parm opcbox.mol2
trajin opcbox.mol2
box x 18.8648434  y 18.4781075  z 19.0062408 alpha 90
check reportfile report.opcbox.dat checkxp
run
EOF
  RunCpptraj "$UNITNAME"
  DoTest report.opcbox.dat.save report.opcbox.dat
fi

EndTest

exit 0
