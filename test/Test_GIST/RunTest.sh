#!/bin/bash

. ../MasterTest.sh

CleanFiles gist.in gist.out gist-*.dx ww_Eij.dat Eww_ij.dat \
           Gist1-*.dx Gist1-*.dat Gist2-*.dx Gist2-*.dat \
           Gist3-*.dx Gist3-*.dat
INPUT="-i gist.in"
TESTNAME='GIST tests'
Requires netcdf notparallel

UNITNAME='Eww Test'
CheckFor notcuda
# doeij test with much smaller grid to save memory
if [ $? -eq 0 ]; then
  cat > gist.in <<EOF
  parm ../tz2.ortho.parm7
  trajin ../tz2.ortho.nc 1 10
  autoimage origin
  gist doorder doeij refdens 0.033422885325 gridcntr 1.44 0.67 0.29 \
    griddim 10 12 10 gridspacn 2.0 prefix Gist1
  go
EOF
  RunCpptraj "GIST water-water interaction test"
  DoTest Gist1-Eww_ij.dat.save Gist1-Eww_ij.dat
fi

# GIST test with finer grid for everything else
UNITNAME='GIST test, orthogonal cell'
cat > gist.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 10
autoimage origin
gist doorder refdens 0.033422885325 gridcntr 1.5 1.0 0.0 \
    griddim 34 44 36 gridspacn 0.50 prefix Gist2 info Info.dat
go
EOF
RunCpptraj "$UNITNAME"
DoTest Gist2-gH.dx.save Gist2-gH.dx
DoTest Gist2-gO.dx.save Gist2-gO.dx
DoTest Gist2-neighbor-norm.dx.save Gist2-neighbor-norm.dx
DoTest Gist2-order-norm.dx.save Gist2-order-norm.dx
# NOTE: gist.out allowed to fail on windows; differences due to slightly
#       difference implementation of printf '%g' (manifests as round-off).
#       THIS IS THE SAVED OUTPUT FROM THE ORIGINAL GIST COMMAND.
DoTest gist.out.save Gist2-output.dat -r 0.0001

# GIST test, nonorthogonal cell
UNITNAME='GIST test, nonorthogonal cell'
cat > gist.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 10
autoimage origin
gist doorder refdens 0.033422885325 \
  gridcntr 0.81 -1.0 0.08 \
  griddim 42 36 40 gridspacn 0.50 \
  prefix Gist3 info Info.dat
EOF
RunCpptraj "$UNITNAME"
DoTest Gist3-gH.dx.save Gist3-gH.dx
DoTest Gist3-gO.dx.save Gist3-gO.dx
DoTest Gist3-neighbor-norm.dx.save Gist3-neighbor-norm.dx
DoTest Gist3-order-norm.dx.save Gist3-order-norm.dx
# NOTE: gist.out allowed to fail on windows; differences due to slightly
#       difference implementation of printf '%g' (manifests as round-off).
DoTest Gist3-output.dat.save Gist3-output.dat -r 0.0001

EndTest
exit 0
