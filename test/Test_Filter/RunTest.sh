#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles filter.in filter.crd filter.dat datafilter.dat

CheckNetcdf

# Test 1
TESTNAME="Filter test."
NotParallel "$TESTNAME"
if [ "$?" -ne 0 ] ; then
  EndTest
  exit 0
fi
TOP='../tz2.truncoct.parm7'
INPUT='filter.in'
cat > filter.in <<EOF
trajin ../tz2.truncoct.nc
rms R1 first :2-11
filter R1 min 0.7 max 0.8 out filter.dat
outtraj filter.crd 
EOF
RunCpptraj "$TESTNAME"
DoTest ../Test_Outtraj/maxmin.crd.save filter.crd
DoTest filter.dat.save filter.dat

# Data filter test
cat > filter.in <<EOF
readdata ../Test_Diffusion/diff_a.xmgr index 1 name A as dat
datafilter A:2 A:3 min 0.0 max 1.2 out datafilter.dat multi name FA
EOF
RunCpptraj "Data filter test."
DoTest datafilter.dat.save datafilter.dat

EndTest

exit 0
