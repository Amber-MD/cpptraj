#!/bin/bash

. ../MasterTest.sh

in=density.in
out1=number_density.dat
out2=mass_density.dat
out3=charge_density.dat
out4=electron_density.dat

CleanFiles $in $out1 $out2 $out3 $out4

NotParallel "Density test."
if [[ $? -ne 0 ]] ; then
  echo ""
  exit 0
fi

INPUT="-i $in"

del='delta 0.25'
masks='":PC@P31" ":PC@N31" ":PC@C2" ":PC | :OL | :OL2"'

cat > $in <<EOF
# crd/top courtesy of Callum Dickson, Imperial College London
parm ../DOPC.parm7
trajin ../DOPC.rst7

center ":PC | :OL | :OL2" origin

density out $out1 number $del $masks
density out $out2 mass $del $masks
density out $out3 charge $del $masks
density out $out4 electron $del $masks
EOF

RunCpptraj "Density Test."
DoTest ${out1}.save $out1
DoTest ${out2}.save $out2
DoTest ${out3}.save $out3
DoTest ${out4}.save $out4
CheckTest
EndTest

exit 0
