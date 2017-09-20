#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in cd.dat cf.dat ci.dat cvt.dat
TESTNAME='clusterdihedral test'
Requires netcdf notparallel
# clusterdihedral
TOP="../tz2.parm7"
INPUT="ptraj.in"
cat > ptraj.in <<EOF
noprogress
trajin ../tz2.nc
clusterdihedral TZ2 :2-12 phibins 2 psibins 2 \
  out cd.dat framefile cf.dat clusterinfo ci.dat clustervtime cvt.dat
EOF
RunCpptraj "$TESTNAME"
DoTest cd.dat.save cd.dat
DoTest cf.dat.save cf.dat
DoTest ci.dat.save ci.dat
DoTest cvt.dat.save cvt.dat

EndTest

exit 0
