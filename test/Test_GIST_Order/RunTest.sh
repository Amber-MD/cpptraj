#!/bin/bash

. ../MasterTest.sh

CleanFiles gist.in gist-*.dx gistout.dat tip4p.dat

TESTNAME='GIST tetrahedral water cluster tests'
Requires notparallel

INPUT="-i gist.in"

UNITNAME='TIP3P water cluster'
cat > gist.in <<EOF
parm test_gist.prmtop
trajin test_gist-center.crd
gist doorder gridcntr 0.5 0.5 0.5 griddim 1 1 1 gridspacn 1.0 out gistout.dat refdens 0.0333
EOF
RunCpptraj "$UNITNAME"
DoTest gist-order-norm.dx.save gist-order-norm.dx
DoTest gistout.dat.save gistout.dat

UNITNAME='TIP4P water cluster'
cat > gist.in <<EOF
parm tip4pew.cluster.parm7
trajin tip4pew.cluster.rst7
gist doorder gridcntr 0.5 0.5 0.5 griddim 1 1 1 gridspacn 1.0 out tip4p.dat refdens 0.0333
EOF
RunCpptraj "$UNITNAME"
#DoTest gist-order-norm.dx.save gist-order-norm.dx
DoTest tip4p.dat.save tip4p.dat

EndTest
exit 0
