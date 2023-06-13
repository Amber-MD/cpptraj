#!/bin/bash

. ../MasterTest.sh

CleanFiles grid.in tyr.rmsfit.dx tyr.rmsfit.dcd tyr.gridfit.dx fitwithreaddata.dx

TESTNAME='Grid rotation tests.'
Requires netcdf

RunTyr() {
  INPUT='-i grid.in'
  cat > grid.in <<EOF
parm dry.tyr.opc3.parm7
noprogress
trajin dry.tyr.nc

set DX = 0.1
set NBINS = 68
set GRID = \$NBINS \$DX \$NBINS \$DX \$NBINS \$DX

# Grid after rms-fitting coordinates
align first :1
grid out tyr.rmsfit.dx :1 name RmsFit \$GRID maskcenter :1
#outtraj tyr.rmsfit.dcd nobox
run

# Fit grid to coordinates. Should only differ from tyr.rmsfit.dx
# in origin location.
grid out tyr.gridfit.dx :1 name GridFit \$GRID rmsfit :1
run
EOF
  RunCpptraj "RMS-fitting grid test"
  DoTest tyr.gridfit.dx.save tyr.gridfit.dx
}

FixWithReadData() {
  INPUT='-i grid.in'
  UNITNAME='RMS-fitting after reading in grid data'
  CheckFor maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > grid.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
readdata ../Test_Grid/out.dx.save name MyGrid
autoimage origin
grid out fitwithreaddata.dx data MyGrid rmsfit :1-13 @CA
run
EOF
    RunCpptraj "$UNITNAME"
    DoTest ../Test_Grid/out.dx.2.save fitwithreaddata.dx
  fi
}

RunTyr
FixWithReadData

EndTest
