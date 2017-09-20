#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in volmap.dx peaks.xyz volmap2.dx volmap3.dx peaks1.xyz

TESTNAME='VMD VolMap Algorithm test'
Requires netcdf maxthreads 10
TOP="../tz2.ortho.parm7"
INPUT="ptraj.in"

cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc
rms first :1-13
center :1-13 mass origin 
volmap volmap3.dx 1.0 1.0 1.0 :WAT@O  \
       radscale 1.36 size 20,20,20 peakcut 0.10 peakfile peaks1.xyz
EOF
RunCpptraj "$TESTNAME"
DoTest volmap3.dx.save volmap3.dx -r 0.00001
DoTest peaks1.xyz.save peaks1.xyz -r 0.00001

UNITNAME='VMD VolMap Algorithm longer tests'
CheckFor long
if [ $? -eq 0 ] ; then
  cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc
rms first :1-13
center :1-13 mass origin 
volmap volmap.dx 0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 \
       radscale 1.36
volmap volmap2.dx 0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 \
       radscale 1.36 peakcut 0.10 peakfile peaks.xyz
EOF
  RunCpptraj "$UNITNAME"
  DoTest volmap.dx.save volmap.dx -r 0.00001
  DoTest volmap.dx.save volmap2.dx -r 0.00001
  DoTest peaks.xyz.save peaks.xyz -r 0.00001
fi

EndTest

exit 0
