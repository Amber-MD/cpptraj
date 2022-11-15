#!/bin/bash

. ../MasterTest.sh

CleanFiles compare.in compare.bonds.dat

TESTNAME='Compare energy tests'
Requires maxthreads 1

INPUT='-i compare.in'

cat > compare.in <<EOF
parm ../tz2.parm7

loadcrd ../tz2.nc 1 1 name CRD1

loadcrd ../tz2.nc 2 2 name CRD2

compareenergy crd0 CRD1 crd1 CRD2 bondout compare.bonds.dat mask1 :2 mask2 :2
EOF
RunCpptraj "$TESTNAME"

EndTest
