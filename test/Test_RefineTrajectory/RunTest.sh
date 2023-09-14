#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in rms.dat

INPUT='-i cpptraj.in'

TESTNAME='Trajectory refinement tests'

Requires netcdf

cat > cpptraj.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
reference ../tz2.nc 1 [FIRST]

# Initial average
rms first !@H=
average crdset R0
run

crdaction R0 rms Rinit ref [FIRST] out rms.dat !@H=
# Loop over references
for i=0;i<10;i++
  j = \$i + 1
  rms ref R\$i !@H=
  average crdset R\$j
  run
  crdaction R\$j rms R\$i.\$j ref R\$i out rms.dat !@H=
done
list
EOF
RunCpptraj "$TESTNAME, using loop"

cat > cpptraj.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc name MyCrd

crdtransform MyCrd mask !@H= 
EOF
RunCpptraj "$TESTNAME, using crdtransform"

EndTest

