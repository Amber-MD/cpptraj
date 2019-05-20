#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in withtime.rst7 withtime.nc

INPUT='-i cpptraj.in'

TESTNAME='Set Time tests'

UNITNAME='Add time to restart test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > cpptraj.in <<EOF
parm ../tz2.parm7
trajin ../tz2.rst7
time time0 10.1
trajout withtime.rst7
go
EOF
  RunCpptraj "$UNITNAME"
  DoTest withtime.rst7.save withtime.rst7
fi

UNITNAME='Add time to trajectory test'
CheckFor netcdf maxthreads 10
if [ $? -eq 0 ] ; then
  cat > cpptraj.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
time time0 100 dt 0.002
trajout withtime.nc
go
EOF
  RunCpptraj "$UNITNAME"
  NcTest withtime.nc.save withtime.nc
fi

EndTest
exit 0
