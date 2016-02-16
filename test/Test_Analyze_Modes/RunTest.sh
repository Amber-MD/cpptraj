#!/bin/bash

. ../MasterTest.sh

CleanFiles modes.in fluct.dat displ.dat corr.dat

TOP=../tz2.parm7
INPUT="modes.in"
CheckPtrajAnalyze
CheckNetcdf

# Test modes fluct and mwcovar matrix generation
TestFluct() {
  cat > modes.in <<EOF
trajin ../tz2.nc
matrix mwcovar name tz2 @CA
analyze matrix tz2 name tz2modes vecs 20
analyze modes fluct stack tz2modes out fluct.dat
EOF
  RunCpptraj "Modes analysis, RMS fluctuations"
  DoTest fluct.dat.save fluct.dat
}

# Test modes displ and modes file read
TestDispl() {
  # Since the displacement test can be fooled by eigenvector
  # sign flips, read in a previously generated set of modes.
  cat > modes.in <<EOF
readdata tz2.evecs.dat name tz2modes
analyze modes displ stack tz2modes out displ.dat
EOF
  RunCpptraj "Modes analysis, displacements"
  DoTest displ.dat.save displ.dat
}

# Test modes corr and mwcovar matrix generation
TestCorr() {
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

TestFluct
TestDispl
TestCorr

EndTest

exit 0
