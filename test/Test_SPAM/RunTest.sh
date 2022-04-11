#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles spam.in spam.dat spam.info spampure.dat test.mdcrd summary.dat

# Check libraries
TESTNAME='SPAM tests'
Requires netcdf maxthreads 10

INPUT="spam.in"
# SPAM test
TOP='../tz2.truncoct.parm7'
if [ -z "$DO_PARALLEL" ] ; then
  cat > spam.in <<EOF
trajin ../tz2.truncoct.nc
autoimage
spam peaks.xyz name SPAM cut 12.0 info spam.info out spam.dat reorder \
     summary summary.dat
trajout test.mdcrd onlyframes 1-2
EOF
  RunCpptraj "SPAM Test"
  DoTest summary.dat.save summary.dat
  DoTest test.mdcrd.save test.mdcrd
  DoTest spam.info.save spam.info
  DoTest spam.dat.save spam.dat
else
  # onlyframes option not supported in parallel
  cat > spam.in <<EOF
trajin ../tz2.truncoct.nc
autoimage
spam peaks.xyz name SPAM cut 12.0 info spam.info out spam.dat reorder \
     summary summary.dat
EOF
  RunCpptraj "SPAM Test"
  DoTest summary.dat.save summary.dat
  DoTest spam.info.save spam.info
  DoTest spam.dat.save spam.dat
fi

# Pure water test
TOP="../spcbox.parm7"
cat > spam.in << EOF
trajin ../spcbox.nc
spam purewater name SPAM cut 12.0 out spampure.dat
EOF
RunCpptraj "SPAM Pure Solvent Test"
DoTest spampure.dat.save spampure.dat

EndTest

exit 0
