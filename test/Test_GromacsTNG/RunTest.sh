#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in temperature.dat rmsd.dat ene.dat

INPUT='-i cpptraj.in'

cat > cpptraj.in <<EOF
parm topol.parm7
trajin md_1_1.tng
# Temperature calc tests velocity read
temperature MyTemp ntc 2 out temperature.dat xmin 0 xstep 10 prec 12.7
# Energy (specifically kinetic VV) tests velocity/force read.
# It is probably not correct since this really requires the half
# step velocities, but is good enough for a regression test.
energy MyEne ^1 out ene.dat xmin 0 xstep 10 prec 12.7 \
  kinetic ketype vv dt 0.002 
# RMSD tests coordinates read
rms MyRms first ^1&@C,CA,N
run
# Convert RMSD to nm
MyRmsNm = MyRms / 10.0
writedata rmsd.dat xmin 0 xstep 10 MyRmsNm prec 12.7
EOF
RunCpptraj "TNG read test"
DoTest temperature.dat.save temperature.dat
DoTest rmsd.dat.save rmsd.dat
DoTest ene.dat.save ene.dat

EndTest

exit 0
