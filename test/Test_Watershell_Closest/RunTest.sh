#!/bin/bash

. ../MasterTest.sh

CleanFiles ws.in wsc.noimage.agr wsc.noimage.avg.dat wsc.closest.dat temp.closest.parm7 temp.closest.nc

TESTNAME='Watershell/Closest tests'
Requires netcdf maxthreads 1
INPUT="ws.in"

# Orthorhombic imaging
#TOP=../tz2.ortho.parm7
#cat > ws.in <<EOF
#trajin ../tz2.ortho.nc
#watershell !:WAT ws.ortho.agr Tz2
#EOF
#RunCpptraj "Watershell Test, orthorhombic imaging."
#DoTest ws.ortho.agr.save ws.ortho.agr

# No imaging
TOP=../tz2.ortho.parm7
cat > ws.in <<EOF
trajin ../tz2.ortho.nc
autoimage
watershell !:WAT out wsc.noimage.agr noimage WS
run
runanalysis avg WS[lower] out wsc.noimage.avg.dat name AvgLower

autoimage
closest \$AvgLower[avg] !:WAT closestout wsc.closest.dat noimage #parmout temp.closest.parm7
#trajout temp.closest.nc
EOF
RunCpptraj "Watershell/Closest Test, no imaging."
DoTest wsc.noimage.agr.save wsc.noimage.agr
DoTest wsc.noimage.avg.dat.save wsc.noimage.avg.dat
DoTest wsc.closest.dat.save wsc.closest.dat

EndTest

exit 0
