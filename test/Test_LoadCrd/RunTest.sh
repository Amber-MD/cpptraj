#!/bin/bash

. ../MasterTest.sh

CleanFiles load.in d1-10.dat d1-12.dat frm.d1-10.dat trj.d1-10.dat \
           trajin.d1-10.dat frm.d1-12.dat

TESTNAME='COORDS data set tests'
Requires netcdf
INPUT="-i load.in"

# Generate comparison file
Generate() {
  cat > load.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
trajin ../tz2.nc 30 88 2
distance d1-10 :1@CA :10@CA out d1-10.dat.save
EOF
  RunCpptraj "Generating comparison file"
}
#Generate

# Loadcrd Test
LoadCrdTest() {
  cat > load.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc CRD
loadcrd ../tz2.nc CRD 30 88 2
crdaction CRD distance d1-10 :1@CA :10@CA out d1-10.dat
list dataset
EOF
  RunCpptraj "LoadCrd test."
  DoTest d1-10.dat.save d1-10.dat
}

# Loadcrd append diff files test
LoadCrdAppend() {
  cat > load.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc TZ2
loadcrd ../tz2.pdb TZ2
loadcrd ../tz2.rst7 TZ2
crdaction TZ2 distance d1-12 :1 :12 out d1-12.dat
list dataset
EOF
  RunCpptraj "LoadCrd append from multiple files test"
  DoTest d1-12.dat.save d1-12.dat
} 

# Loadtraj test
LoadTrajTest() {
  cat > load.in <<EOF
parm ../tz2.parm7
loadtraj ../tz2.nc name CRD
loadtraj ../tz2.nc name CRD 30 88 2
crdaction CRD distance d1-10 :1@CA :10@CA out trj.d1-10.dat
list dataset
EOF
  RunCpptraj "LoadTraj test."
  DoTest d1-10.dat.save trj.d1-10.dat
}

# Loadtraj from trajin test
LoadTrajFromTrajinTest() {
  cat > load.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 
trajin ../tz2.nc 30 88 2
loadtraj name CRD
crdaction CRD distance d1-10 :1@CA :10@CA out trajin.d1-10.dat
list dataset
EOF
  RunCpptraj "LoadTraj from trajin test."
  DoTest d1-10.dat.save trajin.d1-10.dat
}

# Double precision loadcrd Test
LoadFrmTest() {
  cat > load.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc CRD prec double
loadcrd ../tz2.nc CRD 30 88 2 prec double
crdaction CRD distance d1-10 :1@CA :10@CA out frm.d1-10.dat
list dataset
EOF
  RunCpptraj "LoadCrd double precision test."
  DoTest d1-10.dat.save frm.d1-10.dat
}

# Loadcrd append diff files test
LoadFrmAppend() {
  cat > load.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc TZ2 prec double
loadcrd ../tz2.pdb TZ2 prec double
loadcrd ../tz2.rst7 TZ2 prec double
crdaction TZ2 distance d1-12 :1 :12 out frm.d1-12.dat
list dataset
EOF
  RunCpptraj "LoadCrd double precision append from multiple files test"
  DoTest d1-12.dat.save frm.d1-12.dat
}

LoadCrdTest
LoadTrajTest
LoadTrajFromTrajinTest
LoadCrdAppend
LoadFrmTest
LoadFrmAppend

EndTest
exit 0
