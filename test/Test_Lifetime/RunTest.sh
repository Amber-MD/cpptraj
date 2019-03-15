#!/bin/bash

. ../MasterTest.sh

CleanFiles avg.life.5.gnu max.life.5.gnu life.5.gnu life.in perres.avg.gnu \
           perres.cumulative.gnu crv.life.5.gnu solutehb.gnu \
           avg.life2.5.gnu max.life2.5.gnu life2.5.gnu crv.life2.5.gnu \
           solutehb2.gnu
TESTNAME='Lifetime analysis tests'
Requires netcdf

# Basic lifetime test
cat > life.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
hbond H1 series @N,H,C,O #uuseries solutehb.gnu
go
writedata solutehb.gnu \
  H1[solutehb]:24  H1[solutehb]:104 H1[solutehb]:414 H1[solutehb]:422 \
  H1[solutehb]:494 H1[solutehb]:652 H1[solutehb]:184 H1[solutehb]:732 \
  H1[solutehb]:342 H1[solutehb]:352 H1[solutehb]:264 H1[solutehb]:262 \
  H1[solutehb]:424 H1[solutehb]:218 H1[solutehb]:66  H1[solutehb]:188 \
  H1[solutehb]:42  H1[solutehb]:86  H1[solutehb]:794
lifetime H1[solutehb] out life.5.gnu window 5
runanalysis
EOF
INPUT="-i life.in"
RunCpptraj "Lifetime test."
DoTest avg.life.5.gnu.save avg.life.5.gnu
DoTest max.life.5.gnu.save max.life.5.gnu
DoTest life.5.gnu.save life.5.gnu
DoTest crv.life.5.gnu.save crv.life.5.gnu
DoTest solutehb.gnu.save solutehb.gnu

UNITNAME='Lifetime test with data set caching'
CheckFor notparallel
if [ $? -eq 0 ] ; then
  # Using data set caching.
  cat > life.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
usediskcache on
hbond H1 series @N,H,C,O #uuseries solutehb.gnu
go
usediskcache off
writedata solutehb2.gnu title solutehb.gnu \
  H1[solutehb]:24  H1[solutehb]:104 H1[solutehb]:414 H1[solutehb]:422 \
  H1[solutehb]:494 H1[solutehb]:652 H1[solutehb]:184 H1[solutehb]:732 \
  H1[solutehb]:342 H1[solutehb]:352 H1[solutehb]:264 H1[solutehb]:262 \
  H1[solutehb]:424 H1[solutehb]:218 H1[solutehb]:66  H1[solutehb]:188 \
  H1[solutehb]:42  H1[solutehb]:86  H1[solutehb]:794
lifetime H1[solutehb] out life2.5.gnu window 5
datafile avg.life2.5.gnu title avg.life.5.gnu
datafile max.life2.5.gnu title max.life.5.gnu
datafile life2.5.gnu title life.5.gnu
datafile crv.life2.5.gnu title crv.life.5.gnu
EOF
RunCpptraj "$UNITNAME."
  DoTest avg.life.5.gnu.save avg.life2.5.gnu
  DoTest max.life.5.gnu.save max.life2.5.gnu
  DoTest life.5.gnu.save life2.5.gnu
  DoTest crv.life.5.gnu.save crv.life2.5.gnu
  DoTest solutehb.gnu.save solutehb2.gnu
fi

# Average/cumulative average
cat > life.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
rms R1 perres !(@H=) perresmask !(@H=)
go
lifetime R1[res] averageonly out perres.avg.gnu window 5
lifetime R1[res] averageonly cumulative out perres.cumulative.gnu window 5
runanalysis
EOF
RunCpptraj "Lifetime average and cumulative average test."
DoTest perres.avg.gnu.save perres.avg.gnu
DoTest perres.cumulative.gnu.save perres.cumulative.gnu

EndTest

exit 0
