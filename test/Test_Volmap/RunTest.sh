#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in volmap.dx

INPUT="ptraj.in"

# dipole
TOP="../tz2.ortho.parm7"
cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc
rms first :1-13
center :1-13 mass origin 
volmap volmap.dx 0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 \\
       radscale 1.36
volmap volmap2.dx 0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 \\
       radscale 1.36 peakcut 0.10 peakfile peaks.xyz
EOF
RunCpptraj "VMD VolMap Algo. test"
DoTest volmap.dx.save volmap.dx
DoTest volmap.dx.save volmap2.dx
DoTest peaks.xyz.save peaks.xyz
CheckTest

EndTest

exit 0
