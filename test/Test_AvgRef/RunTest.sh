#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rmsavg.in rmsavg.dat 

CheckNetcdf
TOP="../tz2.parm7"

# Test 1
cat > rmsavg.in <<EOF
noprogress
reference ../tz2.nc average  
trajin ../tz2.nc
rms reference :2-11 out rmsavg.dat
EOF
INPUT="rmsavg.in"
RunCpptraj "RMSD Test with averaged reference coordinates."
DoTest rmsavg.dat.save rmsavg.dat

CheckTest

EndTest

exit 0
