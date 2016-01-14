#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles spam.in spam.dat spam.info spampure.dat test.mdcrd

NotParallel "SPAM test."
if [[ $? -ne 0 ]] ; then
  EndTest
  exit 0
fi

# Check libraries
CheckNetcdf

cat > spam.in <<EOF
trajin ../tz2.truncoct.nc
autoimage
spam peaks.xyz name SPAM cut 12.0 info spam.info out spam.dat reorder
trajout test.mdcrd onlyframes 1-2
EOF

INPUT="spam.in"
TOP="../tz2.truncoct.parm7"
RunCpptraj "SPAM Test"

DoTest test.mdcrd.save test.mdcrd
DoTest spam.info.save spam.info
DoTest spam.dat.save spam.dat

cat > spam.in << EOF
trajin ../spcbox.nc
spam purewater name SPAM cut 12.0 out spampure.dat
EOF

INPUT="spam.in"
TOP="../spcbox.parm7"
RunCpptraj "SPAM Pure Solvent Test"

DoTest spampure.dat.save spampure.dat

CheckTest

EndTest

exit 0
