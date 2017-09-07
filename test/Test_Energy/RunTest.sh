#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in ene.agr

INPUT="-i ene.in"

RequiresNetcdf "Simple Energy test"

cat > ene.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
energy DPDP out ene.agr
EOF
RunCpptraj "Simple Energy test"
DoTest ene.agr.save ene.agr

EndTest
exit 0
