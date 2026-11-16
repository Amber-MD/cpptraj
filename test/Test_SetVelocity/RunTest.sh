#!/bin/bash

. ../MasterTest.sh

CleanFiles setvel.in tz2.vel.rst7 V1.dat tz2.scale.rst7 tz2.add.rst7 temp.rst7 \
           t2z.vel.zeromomentum.crd

INPUT='-i setvel.in'

TESTNAME='Set Velocity tests'
Requires maxthreads 1

UNITNAME='Set Velocity test'
cat > setvel.in <<EOF
rng setdefault marsaglia
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
vector V1 momentum magnitude
create V1.dat V1[Mag]
trajout t2z.vel.zeromomentum.crd mdvel
go
EOF
RunCpptraj "$UNITNAME"
DoTest tz2.vel.rst7.save tz2.vel.rst7
DoTest V1.dat.save V1.dat
DoTest t2z.vel.zeromomentum.crd.save t2z.vel.zeromomentum.crd

UNITNAME='Scale velocity test'
cat > setvel.in <<EOF
parm ../tz2.parm7
trajin tz2.vel.rst7.save
strip !:13
setvelocity scale factor 2.0
trajout tz2.scale.rst7
EOF
RunCpptraj "$UNITNAME"
DoTest tz2.scale.rst7.save tz2.scale.rst7

UNITNAME='Add velocity test'
cat > setvel.in <<EOF
parm ../tz2.parm7
trajin tz2.vel.rst7.save
strip !:13
outtraj temp.rst7
setvelocity add value 2.0 vz 0
trajout tz2.add.rst7
EOF
RunCpptraj "$UNITNAME"
DoTest tz2.add.rst7.save tz2.add.rst7

EndTest
exit 0
