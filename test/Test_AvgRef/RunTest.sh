#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rmsavg.in rmsavg.dat tz2.ncrst

TESTNAME='Averaged reference coordinates tests'
Requires netcdf
TOP="../tz2.parm7"

INPUT="rmsavg.in"
if [ -z "$DO_PARALLEL" ] ; then
  cat > rmsavg.in <<EOF
loadcrd ../tz2.nc name TZ2
crdaction TZ2 average crdset Avg_TZ2
trajin ../tz2.nc
rms R0 reference :2-11 out rmsavg.dat
EOF
  RunCpptraj "RMSD Test with averaged reference coordinates."
else
  cat > rmsavg.in <<EOF
trajin ../tz2.nc
average tz2.ncrst
EOF
  RunCpptraj "Parallel, average reference coordinates."
  cat > rmsavg.in <<EOF
trajin ../tz2.nc
reference tz2.ncrst
rms R0 reference :2-11 out rmsavg.dat
EOF
  RunCpptraj "Parallel, RMSD Test with averaged reference coordinates."
fi
DoTest rmsavg.dat.save rmsavg.dat

EndTest

exit 0
