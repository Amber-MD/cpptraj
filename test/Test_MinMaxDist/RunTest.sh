#!/bin/bash

. ../MasterTest.sh

CleanFiles mmdist.in tz2.byatom.dat temp.dat

TESTNAME='Min/max distance tests'

INPUT='-i mmdist.in'
cat > mmdist.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd

mindist name Min mask1 :2 mask2 :11 byatom out tz2.byatom.dat
maxdist name Max mask1 :2 mask2 :11 byatom out tz2.byatom.dat

#nativecontacts :2 :11 savenonnative mindist out temp.dat
#nativecontacts :2 :11 savenonnative maxdist out temp.dat
#precision temp.dat 8 3
EOF
RunCpptraj "Min/max distance test"
DoTest tz2.byatom.dat.save tz2.byatom.dat

EndTest
exit 0
