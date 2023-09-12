#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in rms.dat

INPUT='-i cpptraj.in'

TESTNAME='Trajectory refinement tests'

Requires netcdf

cat > cpptraj.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc

# Initial average
rms first !@H=
average crdset R0
run

# Loop over references
for i=0;i<10;i++
  j = \$i + 1
  rms ref R\$i !@H=
  average crdset R\$j
  run
  crdaction R\$j rms R\$i.\$j ref R\$i out rms.dat
done
list
EOF

RunCpptraj "$TESTNAME"

EndTest

