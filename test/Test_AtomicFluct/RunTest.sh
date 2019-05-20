#!/bin/bash

. ../MasterTest.sh

CleanFiles atomic.in fluct.*.dat dpdp.fluct.dat dpdp.adp.dat \
           fluct.2.pdb occ.2.pdb scale.2.pdb
TESTNAME='Atomic fluctuations tests' 
Requires netcdf

INPUT="atomic.in"

WriteInput() {
  TOP="../tz2.parm7"
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

TOP=../tz2.parm7
cat > $INPUT <<EOF
trajin ../tz2.nc
atomicfluct A0 :2-12
average crdset MyAvg
run
crdout MyAvg fluct.2.pdb bfacdata A0
crdout MyAvg occ.2.pdb occdata A0
crdout MyAvg scale.2.pdb bfacdata A0 bfacscale
EOF
RunCpptraj "Atomicfluct test with PDB B-factor/occupancy output."
DoTest fluct.2.pdb.save fluct.2.pdb
DoTest occ.2.pdb.save occ.2.pdb
DoTest scale.2.pdb.save scale.2.pdb

EndTest

exit 0
