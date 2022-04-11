#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles filter.in filter.crd filter.dat datafilter.dat A2.filtered.dat count.dat

INPUT='filter.in'
TOP='../tz2.truncoct.parm7'

# Test 1
TESTNAME="Filter tests"
Requires notparallel

UNITNAME='Basic filter test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > filter.in <<EOF
trajin ../tz2.truncoct.nc
rms R1 first :2-11
filter R1 min 0.7 max 0.8 out filter.dat
outtraj filter.crd 
EOF
  RunCpptraj "$UNITNAME"
  DoTest ../Test_Outtraj/maxmin.crd.save filter.crd
  DoTest filter.dat.save filter.dat
fi

# Data filter test
cat > filter.in <<EOF
readdata ../Test_Diffusion/diff_a.xmgr.save index 1 name A as dat
datafilter A:2 A:3 min 0.0 max 1.2 out datafilter.dat multi name FA
datafilter A:2 min 0.0 max 1.2 filterset A:2
writedata A2.filtered.dat A:2
list
EOF
RunCpptraj "Data filter test."
DoTest datafilter.dat.save datafilter.dat
DoTest A2.filtered.dat.save A2.filtered.dat

# Test filtering 2D matrix
cat > filter.in <<EOF
readdata ../Test_Matrix/mtest.dat.13.save name Mat read2d square2d
datafilter Mat min .001 max .003 countout count.dat noxcol name Filtered
list
EOF
RunCpptraj "Filter matrix test."
DoTest count.dat.save count.dat

EndTest

exit 0
