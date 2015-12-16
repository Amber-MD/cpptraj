#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles temp.in rmsd.dat
# Test
cat > temp.in <<EOF
noprogress
trajin ../tz2.nc 
rms R0 first out rmsd.dat @%CT
EOF
INPUT="temp.in"
TOP="../tz2.parm7"
RunCpptraj "Selection by atom type test."
DoTest rmsd.dat.save rmsd.dat
CheckTest

EndTest

exit 0
