#!/bin/bash

. ../MasterTest.sh

CleanFiles avg.life.5.gnu max.life.5.gnu life.5.gnu life.in perres.avg.gnu \
           perres.cumulative.gnu crv.life.5.gnu solutehb.gnu

cat > life.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
hbond H1 series @N,H,C,O uuseries solutehb.gnu
go
#create solutehb.gnu H1[solutehb]
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
CheckTest

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
