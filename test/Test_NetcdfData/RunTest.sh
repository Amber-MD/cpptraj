#!/bin/bash

. ../MasterTest.sh

CleanFiles ncdata.in d1.nc d1.dat rmsf.dat ascii.dat.save ascii.dat \
           matrix.nc \
           ca.matrix.nc ca.matrix.dat ca.matrix.dat.save \
           n.ca.matrix.nc n.ca.matrix.dat n.ca.matrix.dat.save \
           ca.rms2d.dat ca.rms2d.dat.save \
           MyEvecs.dat.save MyEvecs.dat WrittenEvecs.dat \
           heavyAtom.matrix.dat.save heavyAtom.matrix.dat \
           heavyEvecs.dat.save heavyEvecs.dat \
           grid.nc grid.dx.save grid.dx \
           hbond.nc hbond.dat.save hbond.dat \
           vector.nc vector.dat.save vector.dat \
           peaks1.nc peaks1.dat.save peaks1.dat \
           cluster.nc C?.matrix.dat.save C?.matrix.dat \
           C1.dat.save C1.dat \
           random.nc random.dat.save random.dat

TESTNAME='NetCDF data file tests.'
Requires netcdf

INPUT='-i ncdata.in'

Write1d() {
  UNITNAME='Write basic 1D NetCDF data'
  SFX='nc'
  CheckFor maxthreads 10
  if [ $? -eq 0 ] ; then
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
    SetConditionalTol d1.nc mpi -a 0.0000000000001
    NcTest d1.nc.save d1.nc $TEST_TOLERANCE
  fi
}

Read1d() {
  UNITNAME='Read basic 1D NetCDF data'
  # Needs 10 procs for save to have been generated in Write1d()
  CheckFor maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > ncdata.in <<EOF
readdata d1.nc.save
list
writedata ascii.dat d1 a1 Dih_00003 Fluct_00004
quit
EOF
    RunCpptraj "$UNITNAME"
    DoTest ascii.dat.save ascii.dat
  fi
}

Write2d() {
  UNITNAME='Write basic 2D NetCDF data'
  cat > ncdata.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
matrix name CA covar @CA out matrix.nc
matrix name N.CA mwcovar @N @CA out matrix.nc
matrix name heavyAtom !@H= mwcovar out matrix.nc
2drms RMS2D @CA out matrix.nc
diagmatrix CA out matrix.nc name MyEvecs
diagmatrix heavyAtom out matrix.nc name heavyEvecs
run
writedata ca.matrix.dat.save CA nosquare2d prec 12.4
writedata n.ca.matrix.dat.save N.CA nosquare2d prec 12.4
writedata heavyAtom.matrix.dat.save heavyAtom nosquare2d prec 12.4
writedata ca.rms2d.dat.save RMS2D nosquare2d prec 12.4
writedata MyEvecs.dat.save MyEvecs
writedata heavyEvecs.dat.save heavyEvecs
EOF
  RunCpptraj "$UNITNAME"

  cat > ncdata.in <<EOF
readdata matrix.nc
writedata ca.matrix.dat CA nosquare2d
writedata n.ca.matrix.dat N.CA nosquare2d
writedata heavyAtom.matrix.dat heavyAtom nosquare2d
writedata ca.rms2d.dat RMS2D nosquare2d
writedata WrittenEvecs.dat MyEvecs
runanalysis diagmatrix CA out MyEvecs.dat name MyEvecs2
runanalysis diagmatrix heavyAtom out heavyEvecs.dat name heavyEvecs2
EOF
  RunCpptraj "Read basic 2D NetCDF data"
  DoTest ca.matrix.dat.save ca.matrix.dat
  DoTest n.ca.matrix.dat.save n.ca.matrix.dat
  DoTest heavyAtom.matrix.dat.save heavyAtom.matrix.dat
  DoTest ca.rms2d.dat.save ca.rms2d.dat
  DoTest MyEvecs.dat.save MyEvecs.dat
  DoTest MyEvecs.dat.save WrittenEvecs.dat
  DoTest heavyEvecs.dat.save heavyEvecs.dat
}

Write3d() {
  UNITNAME='Write basic 3D NetCDF data'
  CheckFor maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > ncdata.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
autoimage origin
bounds ^1 name MyGrid dx 0.5 offset 8
run
list
autoimage origin
grid out grid.nc data MyGrid origin ^1
run
writedata grid.dx.save opendx MyGrid
EOF
    RunCpptraj "$UNITNAME"
    cat > ncdata.in <<EOF
readdata grid.nc name Grid0
writedata grid.dx opendx Grid0
EOF
    RunCpptraj "Read basic 3D NetCDF data"
    DoTest grid.dx.save grid.dx
  fi
}

#Closest() {
#  UNITNAME='Write scalar and string data'
#  cat > ncdata.in <<EOF
#parm ../tz2.ortho.parm7
#trajin ../tz2.ortho.nc
#closest 10 ^1 closestout closest.nc name Closest
#run
#writedata closest.dat.save Closest[*]
#EOF
#  RunCpptraj "$UNITNAME"
#}

Hbond() {
  UNITNAME='Write scalar and string data'
  CheckFor maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > ncdata.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
autoimage origin
hbond HB ^1 solventdonor :WAT solventacceptor :WAT@O out hbond.nc
run
writedata hbond.dat.save HB[*]
EOF
    RunCpptraj "$UNITNAME"
 
    cat > ncdata.in <<EOF
readdata hbond.nc name HB
list
writedata hbond.dat HB[*]
EOF
    RunCpptraj "Read scalar and string data"
    DoTest hbond.dat.save hbond.dat
  fi
}

Vector() {
  UNITNAME='Write vector data'
  CheckFor maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > ncdata.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
vector V1 mask :1 :12 out vector.nc magnitude
run
writedata vector.dat.save V1[*]
EOF
    RunCpptraj "$UNITNAME"

    cat > ncdata.in <<EOF
readdata vector.nc name V1
list
writedata vector.dat V1[*]
EOF
    RunCpptraj "Read vector data"
    DoTest vector.dat.save vector.dat
  fi
}

Volmap() {
  UNITNAME='Write vector/scalar data'
  CheckFor maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > ncdata.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
rms first :1-13
center :1-13 mass origin
volmap name MyMap 1.0 1.0 1.0 :WAT@O \
       radscale 1.36 size 20,20,20 \
       peakcut 0.10 peakfile peaks1.nc
run
writedata peaks1.dat.save MyMap[peaks]
EOF
    RunCpptraj "$UNITNAME"

    cat > ncdata.in <<EOF
readdata peaks1.nc
writedata peaks1.dat MyMap[peaks]
EOF
    RunCpptraj "Read vector/scalar data"
    DoTest peaks1.dat.save peaks1.dat
  fi
}

Cluster() {
  UNITNAME='Write cluster data'
  cat > ncdata.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
#cluster C1 :2-10 clusters 3 epsilon 4.0 \
#  pairdist cluster.nc savepairdist
cluster C1 @CA clusters 5 out cluster.nc \
  sieve 5 bestrep cumulative includesieveincalc \
  pairdist C1matrix pairdistfile cluster.nc savepairdist
cluster C2 :2-10 clusters 3 epsilon 4.0 \
  pairdist C2matrix pairdistfile cluster.nc savepairdist
run

writedata C1.matrix.dat.save C1matrix
writedata C2.matrix.dat.save C2matrix
writedata C1.dat.save C1 
EOF
  RunCpptraj "$UNITNAME"

  cat > ncdata.in <<EOF
readdata cluster.nc
writedata C1.matrix.dat C1matrix
writedata C2.matrix.dat C2matrix
writedata C1.dat C1
list
quit
EOF
  RunCpptraj 'Read cluster data'
  DoTest C1.matrix.dat.save C1.matrix.dat
  DoTest C2.matrix.dat.save C2.matrix.dat
  DoTest C1.dat.save C1.dat
}

Random() {
  UNITNAME='Write unsigned integer data'
  cat > ncdata.in <<EOF
rng setdefault marsaglia createset Marsaglia settype int count 10 seed 10 out random.nc
writedata random.dat.save Marsaglia
EOF
  RunCpptraj "$UNITNAME"

  cat > ncdata.in <<EOF
readdata random.nc
list
writedata random.dat Marsaglia
EOF
  RunCpptraj "Read unsigned integer data"
  DoTest random.dat.save random.dat
}

Write1d
Read1d
Write2d
Write3d
#Closest
Hbond
Vector
Volmap
Cluster
Random

EndTest
