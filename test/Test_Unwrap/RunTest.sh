#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in unwrap.crd unwrap.ortho.crd tz2.ortho.wato.dat tor.dat
TESTNAME='Unwrap tests'
Requires netcdf notparallel

INPUT="ptraj.in"
TOP="../tz2.truncoct.parm7"
cat > ptraj.in <<EOF
trajin ../tz2.truncoct.nc 1 2
unwrap byatom 
trajout unwrap.crd title "Test"
EOF
RunCpptraj "Unwrap non-orthogonal test"
DoTest unwrap.crd.save unwrap.crd

TOP="../tz2.ortho.parm7"
cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc 1 2
unwrap byatom
trajout unwrap.ortho.crd title "Test"
EOF
RunCpptraj "Unwrap orthogonal test"
DoTest unwrap.ortho.crd.save unwrap.ortho.crd

cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc
avgbox MyBox
run

unwrap bymol avgucell MyBox[avg]
diffusion Water :WAT@O out tz2.ortho.wato.dat noimage
run
EOF
RunCpptraj "Unwrap with average box correction test"
DoTest tz2.ortho.wato.dat.save tz2.ortho.wato.dat

cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc
unwrap scheme tor
diffusion TOR :WAT@O out tor.dat noimage
EOF
RunCpptraj "Unwrap with toroidal-view-preserving test"
DoTest ../Test_TorDiff/tor.dat.save tor.dat

EndTest

exit 0
