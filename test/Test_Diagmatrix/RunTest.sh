#!/bin/bash

. ../MasterTest.sh

CleanFiles matrix.in MyThermo.dat

TESTNAME='MW covariance matrix thermo analysis test'
Requires netcdf 

INPUT="-i matrix.in"
cat > matrix.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
matrix MyMatrix mwcovar @1-33&!@H=
diagmatrix MyMatrix thermo outthermo MyThermo.dat
EOF
RunCpptraj "$TESTNAME"
DoTest MyThermo.dat.save MyThermo.dat

EndTest
  
exit 0
