#!/bin/bash

. ../MasterTest.sh

CleanFiles modes.in fluct.dat displ.dat corr.dat

TOP=../tz2.parm7
INPUT="modes.in"
CheckPtrajAnalyze
CheckNetcdf

Test1() {
cat > modes.in <<EOF
trajin ../tz2.nc
matrix mwcovar name tz2 @CA
analyze matrix tz2 name tz2modes vecs 20
analyze modes fluct stack tz2modes out fluct.dat
analyze modes displ stack tz2modes out displ.dat
EOF
RunCpptraj "Analyze Modes Fluct/Displ"
DoTest fluct.dat.save fluct.dat
echo "=================================================================================="
echo "|NOTE: The displacement test is currently disabled due to eigenvector sign flips.|"
echo "=================================================================================="
#DoTest displ.dat.save displ.dat
}

Test2() {
cat > modes.in <<EOF
trajin ../tz2.nc
matrix mwcovar name tz2
analyze matrix tz2 name tz2modes vecs 20
EOF
printf "analyze modes corr stack tz2modes out corr.dat" >> modes.in
for ATOMN in 14 38 52 76 91 105 112 134 158 172 196 218 ; do
  ((ATOMH = ATOMN + 1))
  printf " maskp @$ATOMN @$ATOMH" >> modes.in
done
printf "\n" >> modes.in
RunCpptraj "Analyze Modes Corr"
DoTest corr.dat.save corr.dat
}

Test1
Test2

EndTest

exit 0
