#!/bin/bash

. ../MasterTest.sh

CleanFiles gist.in gist.out gist-*.dx ww_Eij.dat Eww_ij.dat \
           Gist1-*.dx Gist1-*.dat Gist2-*.dx Gist2-*.dat \
           Gist3-*.dx Gist3-*.dat Gist4-*.dx Gist4-*.dat \
           Gist5-*.dx Gist5-*.dat
INPUT="-i gist.in"
TESTNAME='GIST tests'
Requires netcdf notparallel

UNITNAME='Eww Test'
CheckFor notcuda
# doeij test with much smaller grid to save memory
if [ $? -eq 0 ]; then
  # Default test tolerance
  TEST_TOLERANCE='0.0001'
  cat > gist.in <<EOF
  parm ../tz2.ortho.parm7
  trajin ../tz2.ortho.nc 1 10
  autoimage origin
  gist nopme doorder doeij refdens 0.033422885325 gridcntr 1.44 0.67 0.29 \
    griddim 10 12 10 gridspacn 2.0 prefix Gist1
  go
EOF
  RunCpptraj "GIST water-water interaction test"
  DoTest Gist1-Eww_ij.dat.save Gist1-Eww_ij.dat
else
  # GPU test tolerance
  TEST_TOLERANCE='0.0003'
fi

# GIST test with finer grid for everything else
UNITNAME='GIST test, orthogonal cell'
cat > gist.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 10
autoimage origin
gist nopme doorder refdens 0.033422885325 gridcntr 1.5 1.0 0.0 \
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
DoTest gist.out.save Gist2-output.dat -a $TEST_TOLERANCE 

# GIST test, nonorthogonal cell
UNITNAME='GIST test, nonorthogonal cell'
cat > gist.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 10
autoimage origin
gist nopme doorder refdens 0.033422885325 \
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
# The CUDA code requires slightly bigger relative error:
### Maximum absolute error in matching lines = 1.10e-05 at line 5484 field 14
### Maximum relative error in matching lines = 1.08e-04 at line 5484 field 14
### Maximum absolute error in matching lines = 2.00e-04 at line 8015 field 17
### Maximum relative error in matching lines = 1.83e-05 at line 35522 field 17
DoTest Gist3-output.dat.save Gist3-output.dat -a $TEST_TOLERANCE 

UNITNAME='PME-GIST test on orthogonal cell'
CheckFor libpme notcuda
if [ $? -eq 0 ] ; then
  cat > gist.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 10
autoimage origin
#gist pme refdens 0.033422885325 gridcntr 17 20 18 griddim 80 90 80 prefix Gist4 info Info.dat
gist pme refdens 0.033422885325 gridcntr 1.5 1.0 0.0 \
    griddim 34 44 36 gridspacn 0.50 prefix Gist4 info Info.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest Gist4-Solute-Etot-pme-dens.dx.save Gist4-Solute-Etot-pme-dens.dx
  DoTest Gist4-Water-Etot-pme-dens.dx.save Gist4-Water-Etot-pme-dens.dx
  DoTest Gist2-gH.dx.save Gist4-gH.dx
  DoTest Gist2-gO.dx.save Gist4-gO.dx
  DoTest Gist2-neighbor-norm.dx.save Gist4-neighbor-norm.dx
  DoTest Gist4-Esw-dens.dx.save Gist4-Esw-dens.dx
  DoTest Gist4-Eww-dens.dx.save Gist4-Eww-dens.dx
  DoTest Gist4-Info.dat.save Gist4-Info.dat -a $TEST_TOLERANCE
  ## Not including this save on the remote repo bc it is too big.
  #if [ -f 'Gist4-output.dat.save' ] ; then
  #  DoTest Gist4-output.dat.save Gist4-output.dat
  #fi
fi

UNITNAME='PME-GIST test on non-orthogonal cell'
CheckFor libpme notcuda
if [ $? -eq 0 ] ; then
  cat > gist.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 10
autoimage origin
#gist pme refdens 0.033422885325 gridcntr 21 21 21 griddim 90 90 90 prefix Gist5 info Info.dat
gist pme refdens 0.033422885325 gridcntr 0.81 -1.0 0.08 \
  griddim 42 36 40 gridspacn 0.50 prefix Gist5 info Info.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest Gist5-Solute-Etot-pme-dens.dx.save Gist5-Solute-Etot-pme-dens.dx
  DoTest Gist5-Water-Etot-pme-dens.dx.save Gist5-Water-Etot-pme-dens.dx
  DoTest Gist3-gH.dx.save Gist5-gH.dx
  DoTest Gist3-gO.dx.save Gist5-gO.dx
  DoTest Gist3-neighbor-norm.dx.save Gist5-neighbor-norm.dx
  DoTest Gist5-Esw-dens.dx.save Gist5-Esw-dens.dx
  DoTest Gist5-Eww-dens.dx.save Gist5-Eww-dens.dx
  DoTest Gist5-Info.dat.save Gist5-Info.dat -a $TEST_TOLERANCE
  # Not including this save on the remote repo bc it is too big.
  #if [ -f 'Gist5-output.dat.save' ] ; then
  #  DoTest Gist5-output.dat.save Gist5-output.dat
  #fi
fi

EndTest
exit 0
