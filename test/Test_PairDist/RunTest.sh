#!/bin/bash

. ../MasterTest.sh

in=Pr.in
out=Pr.dat

CleanFiles $in $out

INPUT="-i $in"

cat > $in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd

pairdist out $out mask "*" delta 0.1
EOF

RunCpptraj "PairDist Test."
DoTest ${out}.save $out
CheckTest
EndTest

exit 0
