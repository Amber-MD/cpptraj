#!/bin/bash

. ../MasterTest.sh

CleanFiles lj.in ene.dat ene.dat.? switch.dat

INPUT='-i lj.in'
TESTNAME='LJ PME tests.'
Requires libpme maxthreads 1

# Basic test
cat > lj.in <<EOF
parm water_2.parm7
trajin water_2.crd

box x 20 y 20 z 20 alpha 90 beta 90 gamma 90
#debug 10
energy out ene.dat prec 16.8 etype pme cut 8.0 dsumtol 0.0000001 skinnb 0.01 ewcoeff 0.3 ewcoefflj 0.3
EOF
RunCpptraj "LJ PME test."
DoTest ene.dat.save ene.dat

# Kappa sweep
i=0
for kappa in '0.25' '0.35' '0.45' '0.5' ; do
  cat > lj.in <<EOF
parm water_2.parm7
trajin water_2.crd

box x 20 y 20 z 20 alpha 90 beta 90 gamma 90
#debug 10
energy out ene.dat.$i prec 16.8 etype pme cut 8.0 dsumtol 0.0000001 skinnb 0.01 \
       ewcoeff 0.3 ewcoefflj $kappa order 8 nfft 64,64,64
EOF
  RunCpptraj "LJ PME kappa sweep test, $kappa"
  if [ $i -gt 0 ] ; then
    DoTest ene.dat.0 ene.dat.$i
  fi
  ((i++))
done

# Test with switching
cat > lj.in <<EOF
parm water_2.parm7
trajin water_2.crd

box x 20 y 20 z 20 alpha 90 beta 90 gamma 90
energy pmeswitch vdw etype pme cut 8.0 dsumtol 0.0000001 skinnb 0.01 \
       ewcoeff 0.3 ljswidth 2.0 ljpme
energy lrswitch  vdw etype pme cut 8.0 dsumtol 0.0000001 skinnb 0.01 \
       ewcoeff 0.3 ljswidth 2.0
run
writedata switch.dat *[vdw] prec 16.8
EOF
RunCpptraj "LJ with switch function test"
DoTest switch.dat.save switch.dat

EndTest
exit 0
