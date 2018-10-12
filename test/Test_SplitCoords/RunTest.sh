#!/bin/bash

. ../MasterTest.sh

CleanFiles splitcoords.in out.crd

TESTNAME='SplitCoords test'
Requires maxthreads 1

INPUT='-i splitcoords.in'
cat > splitcoords.in <<EOF
parm ../tz2.ortho.parm7
loadcrd ../tz2.ortho.rst7 name CRD

crdaction CRD strip !^2-11
splitcoords CRD name Split
crdout Split out.crd
EOF
RunCpptraj "$TESTNAME"
DoTest out.crd.save out.crd

EndTest
exit 0
