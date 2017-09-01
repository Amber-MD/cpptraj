#!/bin/bash

. ../MasterTest.sh

CleanFiles apm.in apm.dat box.dat areaxy.dat mask.dat

INPUT='-i apm.in'
MaxThreads 1 "Area per molecule test."
if [ $? -eq 0 ] ; then
  cat > apm.in <<EOF
parm ../DOPC.parm7
trajin ../DOPC.rst7
areapermol DOPC1 nmols 36 out apm.dat noheader
areapermol DOPC2 :OL nlayers 2 out mask.dat noheader
vector box out box.dat prec 16.8
run

readdata box.dat index 1 name box
areaxy = box:2 * box:3
areaxy = areaxy / 36
writedata areaxy.dat areaxy noheader
quit
EOF
  RunCpptraj "Area per molecule test."
  DoTest apm.dat.save apm.dat
  DoTest apm.dat.save areaxy.dat
  DoTest apm.dat.save mask.dat
fi

EndTest
exit 0
