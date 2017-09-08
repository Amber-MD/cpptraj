#!/bin/bash

. ../MasterTest.sh

CleanFiles systemVF.in test.systemVF.nc test.systemVf.ncdump
TESTNAME='Amber NetCDF trajectory with Velocity/Force test'
Requires netcdf pnetcdf

cat > systemVF.in <<EOF
parm systemVF.parm7
trajin systemVF.nc
trajout test.systemVF.nc netcdf 
EOF
INPUT="-i systemVF.in"
RunCpptraj "$TESTNAME"
NcTest systemVF.nc test.systemVF.nc
EndTest
exit 0
