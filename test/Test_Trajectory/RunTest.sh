#!/bin/bash

. ../MasterTest.sh

TESTNAME='Trajectory format tests'

CleanFiles traj.in

INPUT='-i traj.in'

cat > traj.in <<EOF
parm prmtop/ti.prmtop
parminfo
debug 10
trajin trajin/inpcrd2
trajout test1.rst7
EOF
RunCpptraj "Read/write Amber restart format"

EndTest
