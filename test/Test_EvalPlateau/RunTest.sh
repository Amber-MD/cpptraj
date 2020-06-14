#!/bin/bash

. ../MasterTest.sh

TESTNAME='Evalulate plateau test'

INPUT='-i evalp.in'

CleanFiles evalp.in Eval.agr Eval.results

cat > evalp.in <<EOF
readdata density.dat index 1 name MD
runanalysis evalplateau MD name EQ out Eval.agr resultsout Eval.results
EOF
RunCpptraj "$TESTNAME"
DoTest Eval.results.save Eval.results
DoTest Eval.agr.save Eval.agr

EndTest
