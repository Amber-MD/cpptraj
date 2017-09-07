#!/bin/bash

. ../MasterTest.sh

CleanFiles drmsd.in drmsd.dat

RequiresNetcdf "Distance RMSD test"
cat > drmsd.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
drmsd drms_nofit out drmsd.dat
rms rms_nofit out drmsd.dat nofit
rms rms_fit out drmsd.dat
drmsd drms_fit out drmsd.dat
EOF
INPUT="-i drmsd.in"

RunCpptraj "Distance RMSD test."

DoTest drmsd.dat.save drmsd.dat
EndTest

exit 0
