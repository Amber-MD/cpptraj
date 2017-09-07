#!/bin/bash

. ../MasterTest.sh

CleanFiles gist.in gist-*.dx gistout.dat

RequiresNotParallel "GIST tetrahedral water cluster test."

INPUT="-i gist.in"

cat > gist.in <<EOF
parm test_gist.prmtop
trajin test_gist-center.crd
gist doorder gridcntr 0.5 0.5 0.5 griddim 1 1 1 gridspacn 1.0 out gistout.dat refdens 0.0333
EOF
RunCpptraj "GIST tetrahedral water cluster test."

DoTest gist-order-norm.dx.save gist-order-norm.dx
DoTest gistout.dat.save gistout.dat

EndTest
exit 0
