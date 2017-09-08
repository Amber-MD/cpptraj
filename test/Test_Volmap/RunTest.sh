#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in volmap.dx peaks.xyz volmap2.dx volmap3.dx

TESTNAME='VMD VolMap Algorithm test'
Requires netcdf maxthreads 10
TOP="../tz2.ortho.parm7"
INPUT="ptraj.in"
cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc
rms first :1-13
center :1-13 mass origin 
volmap volmap.dx 0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 \
       radscale 1.36
volmap volmap2.dx 0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 \
       radscale 1.36 peakcut 0.10 peakfile peaks.xyz
volmap volmap3.dx 1.0 1.0 1.0 :WAT@O  centermask !:1-13 \
       radscale 1.36 size 20,20,20
EOF
RunCpptraj "$TESTNAME"
DoTest volmap.dx.save volmap.dx -r 0.00001
DoTest volmap.dx.save volmap2.dx -r 0.00001
DoTest peaks.xyz.save peaks.xyz -r 0.00001
DoTest volmap3.dx.save volmap3.dx -r 0.00001

EndTest

exit 0
