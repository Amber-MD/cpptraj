#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in ene.agr

INPUT="-i ene.in"

TESTNAME='Simple energy test'
Requires netcdf

cat > ene.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
energy DPDP out ene.agr
EOF
RunCpptraj "$TESTNAME"
DoTest ene.agr.save ene.agr

EndTest
exit 0
