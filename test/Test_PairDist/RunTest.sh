#!/bin/bash

. ../MasterTest.sh

in=Pr.in
out=Pr.dat

CleanFiles $in $out ortho.dat

INPUT="-i $in"

cat > $in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd

pairdist P0 out $out mask "*" delta 0.1 maxdist 36.25
EOF

RunCpptraj "PairDist Test."
DoTest ${out}.save $out

cat > $in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
pairdist P1 out ortho.dat mask :WAT@O delta 0.1 maxdist 20.0
EOF
RunCpptraj "Pairdist test, orthogonal imaging."
DoTest ortho.dat ortho.dat.save

EndTest

exit 0
