#!/bin/bash

. ../MasterTest.sh

CleanFiles ncdata.in d1.nc d1.dat rmsf.dat ascii.dat.save ascii.dat \
           ca.matrix.nc ca.matrix.dat

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
clear trajin
trajin ../tz2.nc
align :2-12&!@H=
rmsf :2,4,9,11&!@H= byres out d1.$SFX
run
writedata ascii.dat.save d1 a1 Dih_00003 Fluct_00004
EOF
RunCpptraj "$UNITNAME"
NcTest d1.nc.save d1.nc

UNITNAME='Read basic 1D NetCDF data'
cat > ncdata.in <<EOF
readdata d1.nc.save
list
writedata ascii.dat d1 a1 Dih_00003 Fluct_00004
quit
EOF
RunCpptraj "$UNITNAME"
DoTest ascii.dat.save ascii.dat

UNITNAME='Write basic 2D NetCDF data'
cat > ncdata.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
matrix name CA covar @CA out ca.matrix.nc
run
writedata ca.matrix.dat CA nosquare2d
EOF
RunCpptraj "$UNITNAME"

EndTest
