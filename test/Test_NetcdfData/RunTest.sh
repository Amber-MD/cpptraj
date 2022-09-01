#!/bin/bash

. ../MasterTest.sh

CleanFiles ncdata.in d1.nc

TESTNAME='NetCDF data file tests.'
Requires netcdf

INPUT='-i ncdata.in'

UNITNAME='Write basic 1D NetCDF data'
cat > ncdata.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
distance d1 :1 :12 out d1.nc
EOF
RunCpptraj "$UNITNAME"

EndTest
