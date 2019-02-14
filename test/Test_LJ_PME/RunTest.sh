#!/bin/bash

. ../MasterTest.sh

CleanFiles lj.in ene.dat

INPUT='-i lj.in'
TESTNAME='LJ PME tests.'
Requires libpme maxthreads 1

cat > lj.in <<EOF
parm water_2.parm7
trajin water_2.crd

box x 20 y 20 z 20 alpha 90 beta 90 gamma 90
debug 10
energy out ene.dat prec 16.8 etype pme cut 8.0 dsumtol 0.0000001 skinnb 0.01 ewcoeff 0.3 ewcoefflj 0.3
EOF
RunCpptraj "LJ PME test."
DoTest ene.dat.save ene.dat

EndTest
exit 0
