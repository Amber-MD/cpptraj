#!/bin/bash

. ../MasterTest.sh

CleanFiles setvel.in tz2.vel.rst7 V1.dat

INPUT='-i setvel.in'

TESTNAME="Set Velocity test."
MaxThreads 1 "$TESTNAME"
if [ "$?" -eq 0 ] ; then
  cat > setvel.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
setvelocity tempi 298.0 ig 10
temperature MyTemp
vector V0 momentum
trajout tz2.vel.rst7
go
printdata MyTemp
printdata V0

clear trajin
trajin tz2.vel.rst7
setvelocity modify zeromomentum
vector V1 momentum out V1.dat
go
EOF
  RunCpptraj "$TESTNAME"
  DoTest tz2.vel.rst7.save tz2.vel.rst7
  DoTest V1.dat.save V1.dat
fi
EndTest
exit 0
