#!/bin/bash

. ../MasterTest.sh

CleanFiles order.in ol?.dat ol3.agr

NotParallel "Lipid Order Parameter Test."
if [[ $? -ne 0 ]] ; then
  EndTest
  exit 0
fi

INPUT="-i order.in"

cat > order.in <<EOF
# crd/top courtesy of Callum Dickson, Imperial College London
parm ../DOPC.parm7
trajin ../DOPC.rst7

lipidorder2 :OL,PC  out ol1.dat
#debug actions 3
lipidorder2 :OL2,PC out ol2.dat
lipidorder2 out ol3.agr xydy
precision ol1.dat 10 7
precision ol2.dat 10 7
precision ol3.agr 10 7
datafile ol1.dat xprec 3.0
datafile ol2.dat xprec 3.0
datafile ol3.agr xprec 3.0
EOF

RunCpptraj "Lipid Order Parameter Test 2."
DoTest ol1.dat.save ol1.dat
DoTest ol2.dat.save ol2.dat
DoTest ol3.agr.save ol3.agr

EndTest

exit 0
