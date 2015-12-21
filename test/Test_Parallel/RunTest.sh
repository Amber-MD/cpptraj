#!/bin/bash

# Test some cpptraj actions in parallel

. ../MasterTest.sh

CleanFiles para.in d1-12.dat a1-6-12.dat phi2.dat rmsd.dat avg.rst7 offset.dat test.nc

INPUT="-i para.in"

Test1() {
cat > para.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
distance D1-12 :1 :12 out d1-12.dat
angle A1-6-12 :1 :6 :12 out a1-6-12.dat
dihedral Phi2 :1@C :2@N :2@CA :2@C out phi2.dat
rms first @CA out rmsd.dat
average avg.rst7
EOF
RunCpptraj "Parallel tests."
DoTest d1-12.dat.save d1-12.dat
DoTest a1-6-12.dat.save a1-6-12.dat
DoTest phi2.dat.save phi2.dat
DoTest rmsd.dat.save rmsd.dat
DoTest avg.rst7.save avg.rst7
}

Test2() {
  CheckPnetcdf "Parallel test with trajectory offset"
  if [[ $? -eq 0 ]] ; then
    cat > para.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc 3 90 4
distance D1-12 :1 :12 out offset.dat
trajout test.nc
EOF
    RunCpptraj "Parallel test with trajectory offset."
    DoTest offset.dat.save offset.dat
    NcTest test.nc.save test.nc
  fi
}

Test1
Test2

EndTest
exit 0
