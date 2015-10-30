#!/bin/bash

. ../MasterTest.sh

CleanFiles systemVF.in test.systemVF.nc test.systemVf.ncdump
CheckNetcdf

cat > systemVF.in <<EOF
parm systemVF.parm7
trajin systemVF.nc
trajout test.systemVF.nc netcdf 
EOF
INPUT="-i systemVF.in"
RunCpptraj "Amber Netcdf with Velocity/Force test"
NcTest systemVF.nc test.systemVF.nc
CheckTest
EndTest
exit 0
