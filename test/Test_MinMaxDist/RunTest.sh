#!/bin/bash

. ../MasterTest.sh

CleanFiles mmdist.in tz2.byatom.dat temp.dat tz2.byatom.1.dat \
           tz2.byres.?.dat tz2.bymol.?.dat

TESTNAME='Min/max distance tests'

INPUT='-i mmdist.in'
cat > mmdist.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd

mindist name Min mask1 :2 mask2 :11 byatom out tz2.byatom.dat
maxdist name Max mask1 :2 mask2 :11 byatom out tz2.byatom.dat
mindist name Min1 mask1 ^1 byatom out tz2.byatom.1.dat
maxdist name Max1 mask1 ^1 byatom out tz2.byatom.1.dat


#nativecontacts :2 :11 savenonnative mindist out temp.dat
#nativecontacts :2 :11 savenonnative maxdist out temp.dat
#nativecontacts ^1 maxdist out temp.dat
#precision temp.dat 8 3
EOF
RunCpptraj "Min/max distance by atom tests"
DoTest tz2.byatom.dat.save tz2.byatom.dat
DoTest tz2.byatom.1.dat.save tz2.byatom.1.dat

cat > mmdist.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.crd

mindist name ResMin1 mask1 :1-7 byres out tz2.byres.1.dat
maxdist name ResMax1 mask1 :1-7 byres out tz2.byres.1.dat
mindist name ResMin2 mask1 :2-4 mask2 :9-11 byres out tz2.byres.2.dat
maxdist name ResMax2 mask1 :2-4 mask2 :9-11 byres out tz2.byres.2.dat
EOF
RunCpptraj "Min/max distance by residue tests"
DoTest tz2.byres.1.dat.save tz2.byres.1.dat
DoTest tz2.byres.2.dat.save tz2.byres.2.dat

cat > mmdist.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.crd

mindist name MolMin1 mask1 ^1-10 bymol out tz2.bymol.1.dat
maxdist name MolMax1 mask1 ^1-10 bymol out tz2.bymol.1.dat
mindist name MolMin2 mask1 ^1 mask2 ^342,567,638,1026,1560,1846&@O bymol out tz2.bymol.2.dat
maxdist name MolMax2 mask1 ^1 mask2 ^342,567,638,1026,1560,1846&@O bymol out tz2.bymol.2.dat
#closest 1 ^1 first closestout temp.dat
EOF
RunCpptraj "Min/max distance by molecule tests"
DoTest tz2.bymol.1.dat.save tz2.bymol.1.dat
DoTest tz2.bymol.2.dat.save tz2.bymol.2.dat
EndTest
exit 0
