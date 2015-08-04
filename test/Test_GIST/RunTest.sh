#!/bin/bash

. ../MasterTest.sh

CleanFiles gist.in gist.out gist-*.dx ww_Eij.dat Eww_ij.dat

INPUT="-i gist.in"

cat > gist.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 10
autoimage origin
gist doorder doeij refdens 0.033422885325 gridcntr 1.5 1.0 0.0 \
  griddim 34 44 36 gridspacn 0.50 out gist.out
go
EOF
RunCpptraj "GIST test"
DoTest Eww_ij.dat.save Eww_ij.dat
DoTest gist-gH.dx.save gist-gH.dx
DoTest gist-gO.dx.save gist-gO.dx
DoTest gist-neighbor-norm.dx.save gist-neighbor-norm.dx
DoTest gist-order-norm.dx.save gist-order-norm.dx
DoTest gist.out.save gist.out
EndTest
exit 0
