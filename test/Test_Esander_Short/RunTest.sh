#!/bin/bash

. ../MasterTest.sh

CleanFiles esander.in gb.dat pme.dat

TESTNAME='Energy via sander API short tests (GB/PME)'
Requires sanderlib maxthreads 1

INPUT='-i esander.in'
cat > esander.in <<EOF
parm ../Test_Hbond/strip.4lztSc_nowat.parm7
trajin ../Test_Hbond/strip.4lztSc.rst7
esander 4lztSc out gb.dat \
  saltcon 1.0 cut 9999.0 gbsa 1 igb 1 ntb 0
esander PME out pme.dat
EOF
RunCpptraj "$TESTNAME"
DoTest gb.dat.save gb.dat
DoTest pme.dat.save pme.dat

EndTest
exit 0
