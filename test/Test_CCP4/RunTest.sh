#!/bin/bash

. ../MasterTest.sh

CleanFiles data.in fav8.dx

INPUT="-i data.in"

cat > data.in <<EOF
readdata fav8.guv.O.1.ccp4 name MyCcp4
list dataset
writedata fav8.dx MyCcp4
quit
EOF
RunCpptraj "CCP4 -> DX conversion test."
DoTest fav8.dx.save fav8.dx

EndTest
exit 0
