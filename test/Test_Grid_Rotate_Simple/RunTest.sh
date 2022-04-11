#!/bin/bash

. ../MasterTest.sh

CleanFiles grid.in tyr.rmsfit.dx tyr.rmsfit.dcd tyr.gridfit.dx

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

RunTyr

EndTest
