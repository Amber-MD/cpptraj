#!/bin/bash

. ../MasterTest.sh

TESTNAME='Test support for Lennard-Jones C coefficients'

CleanFiles ljc.in

INPUT='-i ljc.in'

cat > ljc.in <<EOF
parm znf_1264.prmtop.save
parmwrite out cpptraj.znf_1264.parm7
quit
EOF
RunCpptraj "Test writing topology with LJ C coefficients"

EndTest
