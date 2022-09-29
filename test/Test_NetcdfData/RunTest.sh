#!/bin/bash

. ../MasterTest.sh

CleanFiles ncdata.in d1.nc d1.dat

TESTNAME='NetCDF data file tests.'
Requires netcdf

INPUT='-i ncdata.in'

UNITNAME='Write basic 1D NetCDF data'
SFX='nc'
cat > ncdata.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
distance d1 :1 :12 out d1.$SFX
angle a1 :1 :2 :3 out d1.$SFX
run
trajin ../tz2.nc 11 15
dihedral :1 :2 :3 :4 out d1.$SFX
run
EOF
RunCpptraj "$UNITNAME"
NcTest d1.nc.save d1.nc

EndTest
