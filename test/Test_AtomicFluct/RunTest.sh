#!/bin/bash

. ../MasterTest.sh

CleanFiles atomic.in fluct.*.dat 
CheckNetcdf
INPUT="atomic.in"
TOP="../tz2.parm7"

WriteInput() {
  cat > $INPUT <<EOF
trajin ../tz2.nc
atomicfluct out fluct.$2.dat $1
EOF
  RunCpptraj "Cpptraj: Atomic fluctuations test [$1]."
  DoTest fluct.$2.dat.save fluct.$2.dat
  CheckTest 
}

WriteInput "@C,CA,N byres bfactor" 1
WriteInput ":2-12 byatom" 2
WriteInput "bymask :3,4,5" 3
WriteInput "start 10 stop 30 offset 2 byres bfactor" 4

EndTest

exit 0
