#!/bin/bash

. ../MasterTest.sh

CleanFiles ncdata.in d1.nc d1.dat rmsf.dat ascii.dat.save ascii.dat \
           matrix.nc \
           ca.matrix.nc ca.matrix.dat \
           n.ca.matrix.nc n.ca.matrix.dat \
           ca.rms2d.dat \
           MyEvecs.dat.save MyEvecs.dat

TESTNAME='NetCDF data file tests.'
Requires netcdf

INPUT='-i ncdata.in'

Write1d() {
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
}

Read1d() {
  UNITNAME='Read basic 1D NetCDF data'
  cat > ncdata.in <<EOF
readdata d1.nc.save
list
writedata ascii.dat d1 a1 Dih_00003 Fluct_00004
quit
EOF
  RunCpptraj "$UNITNAME"
  DoTest ascii.dat.save ascii.dat
}

Write2d() {
  UNITNAME='Write basic 2D NetCDF data'
  cat > ncdata.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
matrix name CA covar @CA out matrix.nc
matrix name N.CA mwcovar @N @CA out matrix.nc
2drms RMS2D @CA out matrix.nc
diagmatrix CA out MyEvecs.dat.save name MyEvecs
run
writedata ca.matrix.dat.save CA nosquare2d prec 12.4
writedata n.ca.matrix.dat.save N.CA nosquare2d prec 12.4
writedata ca.rms2d.dat.save RMS2D nosquare2d prec 12.4
EOF
  RunCpptraj "$UNITNAME"

  cat > ncdata.in <<EOF
readdata matrix.nc
writedata ca.matrix.dat CA nosquare2d
writedata n.ca.matrix.dat N.CA nosquare2d
writedata ca.rms2d.dat RMS2D nosquare2d
runanalysis diagmatrix CA out MyEvecs.dat name MyEvecs
EOF
  RunCpptraj "Read basic 2D NetCDF data"
  DoTest ca.matrix.dat.save ca.matrix.dat
  DoTest n.ca.matrix.dat.save n.ca.matrix.dat
  DoTest ca.rms2d.dat.save ca.rms2d.dat
  DoTest MyEvecs.dat.save MyEvecs.dat
}

Write1d
Read1d
Write2d

EndTest
