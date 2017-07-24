#!/bin/bash

. ../MasterTest.sh

CleanFiles order.in 

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

lipidorder2 :OL,PC
lipidorder2 :OL2,PC
EOF

RunCpptraj "Lipid Order Parameter Test 2."

EndTest

exit 0
