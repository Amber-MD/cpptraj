#!/bin/bash

. ../MasterTest.sh

CleanFiles apm.in apm.dat box.dat areaxy.dat mask.dat mask2.dat \
           xz.dat areaxz.dat yz.dat areayz.dat

INPUT='-i apm.in'
TESTNAME='Area per molecule test'
Requires maxthreads 1
cat > apm.in <<EOF
parm ../DOPC.parm7
trajin ../DOPC.rst7
areapermol DOPC1 nmols 36 out apm.dat noheader
areapermol DOPC2 :OL nlayers 2 out mask.dat noheader
areapermol DOPC3 :OL,PC,OL2 nlayers 2 out mask2.dat noheader
areapermol DOPC4 nmols 36 xz out xz.dat noheader
areapermol DOPC5 nmols 36 yz out yz.dat noheader
vector box out box.dat prec 16.8
run

readdata box.dat index 1 name box
areaxy = box:2 * box:3
areaxy = areaxy / 36
writedata areaxy.dat areaxy noheader
areaxz = box:2 * box:4
areaxz = areaxz / 36
writedata areaxz.dat areaxz noheader
areayz = box:3 * box:4
areayz = areayz / 36
writedata areayz.dat areayz noheader
quit
EOF
RunCpptraj "$TESTNAME"
DoTest apm.dat.save apm.dat
DoTest apm.dat.save areaxy.dat
DoTest apm.dat.save mask.dat
DoTest apm.dat.save mask2.dat
DoTest areaxz.dat xz.dat
DoTest areayz.dat yz.dat

EndTest
exit 0
