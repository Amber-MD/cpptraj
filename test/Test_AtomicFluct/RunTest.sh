#!/bin/bash

. ../MasterTest.sh

CleanFiles atomic.in fluct.*.dat dpdp.fluct.dat dpdp.adp.dat
TESTNAME='Atomic fluctuations tests' 
Requires netcdf
INPUT="atomic.in"
TOP="../tz2.parm7"

WriteInput() {
  cat > $INPUT <<EOF
trajin ../tz2.nc
atomicfluct out fluct.$2.dat $1
EOF
  RunCpptraj "Atomic fluctuations test [$1]."
  DoTest fluct.$2.dat.save fluct.$2.dat
}

WriteInput "@C,CA,N byres bfactor" 1
WriteInput ":2-12 byatom" 2
WriteInput "bymask :3,4,5" 3
WriteInput "start 10 stop 30 offset 2 byres bfactor" 4

TOP=../DPDP.parm7
cat > $INPUT <<EOF
trajin ../DPDP.nc
rms first mass
atomicfluct out dpdp.fluct.dat adpout dpdp.adp.dat
EOF
RunCpptraj "Atomicfluct test with ADP output"
DoTest dpdp.fluct.dat.save dpdp.fluct.dat
DoTest dpdp.adp.dat.save dpdp.adp.dat

EndTest

exit 0
