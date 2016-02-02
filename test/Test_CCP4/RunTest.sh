#!/bin/bash

. ../MasterTest.sh

CleanFiles data.in fav8.dx fav8.ccp4

INPUT="-i data.in"

cat > data.in <<EOF
readdata fav8.guv.O.1.ccp4 name MyCcp4
list dataset
writedata fav8.dx MyCcp4
writedata fav8.ccp4 MyCcp4 title "Amber 3D-RISM CCP4 map volumetric data. Format revision A."
quit
EOF
RunCpptraj "CCP4 read/write and DX conversion test."
DoTest fav8.dx.save fav8.dx
DoTest fav8.guv.O.1.ccp4 fav8.ccp4

EndTest
exit 0
