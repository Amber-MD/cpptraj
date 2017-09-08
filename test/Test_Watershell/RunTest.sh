#!/bin/bash

. ../MasterTest.sh

CleanFiles ws.in ws.agr ws.ortho.agr ws.noimage.agr

TESTNAME='Watershell tests'
Requires netcdf
INPUT="ws.in"

# Non-orthorhombic imaging
TOP=../tz2.truncoct.parm7
cat > ws.in <<EOF
trajin ../tz2.truncoct.nc
watershell !:WAT ws.agr Tz2
EOF
RunCpptraj "Watershell Test, non-orthorhombic imaging."
DoTest ws.agr.save ws.agr

# Orthorhombic imaging
TOP=../tz2.ortho.parm7
cat > ws.in <<EOF
trajin ../tz2.ortho.nc
watershell !:WAT ws.ortho.agr Tz2
EOF
RunCpptraj "Watershell Test, orthorhombic imaging."
DoTest ws.ortho.agr.save ws.ortho.agr

# No imaging
TOP=../tz2.truncoct.parm7
cat > ws.in <<EOF
trajin ../tz2.truncoct.nc
watershell !:WAT ws.noimage.agr Tz2 noimage
EOF
RunCpptraj "Watershell Test, no imaging."
DoTest ws.noimage.agr.save ws.noimage.agr

EndTest

exit 0
