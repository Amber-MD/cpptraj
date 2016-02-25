#!/bin/bash

. ../MasterTest.sh

in=Pr.in
out=Pr.dat

CleanFiles $in $out

INPUT="-i $in"

cat > $in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd

pairdist P0 out $out mask "*" delta 0.1
EOF

RunCpptraj "PairDist Test."
if [[ -z $DO_PARALLEL ]] ; then
  DoTest ${out}.save $out
else
  # NOTE: In parallel the differences can be larger than expected due
  #       to the way the pairdist command accumulates the mean and
  #       standard deviation.
  DoTest ${out}.save $out -r 0.99
fi
EndTest

exit 0
