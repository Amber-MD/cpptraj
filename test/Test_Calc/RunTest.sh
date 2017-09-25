#!/bin/bash

. ../MasterTest.sh

CleanFiles calc.in calc.dat mag.dat

INPUT='-i calc.in'
cat > calc.in <<EOF
     V0 = 7--sqrt(4)
calc V1 = (5.031661 + 4.334282 + 3.947999) / 3
calc V2 = 5-3+7 prec 3.1
calc V3 = 5-3-7 format general
calc V4 = ln( 10 ) prec 16.8 format scientific
     V5 = (20 - 7) + 3 * (21 / 3)
     V6 = (20 - 7) + 3 * (21 / 3)^2
writedata calc.dat xprec 2.0 invert V*
EOF
RunCpptraj "Calc test"
DoTest calc.dat.save calc.dat

cat > calc.in <<EOF
readdata ../Test_Rotdif/rvecs.dat.save name Vec index 1
X2     = Vec:2 * Vec:2
Y2pZ2  =  (Vec:3^2) + (Vec:4 * Vec:4)
Mag    = sqrt(X2 + Y2pZ2)
AvgMag = avg(Mag)
writedata mag.dat AvgMag
EOF
RunCpptraj "Calc test with data sets"
DoTest mag.dat.save mag.dat

EndTest
exit 0
