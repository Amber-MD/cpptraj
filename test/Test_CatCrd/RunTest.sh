#!/bin/bash

. ../MasterTest.sh

CleanFiles cat.in cat.crd.save cat.crd

TESTNAME='Concatenate COORDS data set test'
Requires maxthreads 1

INPUT='-i cat.in'

cat > cat.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 4
trajout cat.crd.save
go

loadcrd ../tz2.nc 1 2 name frame12
loadcrd ../tz2.nc 3 4 name frame34
catcrd frame12 frame34 name frame1234
crdout frame1234 cat.crd 
quit
EOF
RunCpptraj "$TESTNAME"
DoTest cat.crd.save cat.crd

EndTest
exit 0
