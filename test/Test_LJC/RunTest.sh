#!/bin/bash

. ../MasterTest.sh

TESTNAME='Test support for Lennard-Jones C coefficients'

CleanFiles ljc.in cpptraj.znf_1264.parm7

INPUT='-i ljc.in'

cat > ljc.in <<EOF
parm znf_1264.prmtop.save
parminfo
parmwrite out cpptraj.znf_1264.parm7
quit
EOF
RunCpptraj "Test reading/writing topology with LJ C coefficients"
DoTest cpptraj.znf_1264.parm7.save cpptraj.znf_1264.parm7 -I %VERSION

EndTest
