#!/bin/bash

. ../MasterTest.sh

CleanFiles setvel.in tz2.vel.rst7

INPUT='-i setvel.in'

TESTNAME="Set Velocity test."
MaxThreads 1 "$TESTNAME"
if [ "$?" -eq 0 ] ; then
  cat > setvel.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
setvelocity tempi 298.0 ig 10
temperature MyTemp
trajout tz2.vel.rst7
go
printdata MyTemp
EOF
  RunCpptraj "$TESTNAME"
  DoTest tz2.vel.rst7.save tz2.vel.rst7
fi
EndTest
exit 0
