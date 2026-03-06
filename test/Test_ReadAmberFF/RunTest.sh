#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in tip3pf.parm.dat toyrna.parm.dat toyrna.tip3pf.parm.dat

INPUT='-i cpptraj.in'

cat > cpptraj.in <<EOF
readdata frcmod.tip3pf as frcmod name TIP3PF
writedata tip3pf.parm.dat TIP3PF

readdata toyrna.dat as amberff name TOYRNA
writedata toyrna.parm.dat TOYRNA

readdata toyrna.dat as amberff name PARM
readdata frcmod.tip3pf as frcmod name PARM
writedata toyrna.tip3pf.parm.dat PARM
EOF
RunCpptraj "Test read Amber FF and force modification files."
DoTest tip3pf.parm.dat.save tip3pf.parm.dat
DoTest toyrna.parm.dat.save toyrna.parm.dat
DoTest toyrna.tip3pf.parm.dat.save toyrna.tip3pf.parm.dat

EndTest
