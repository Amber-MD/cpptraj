#!/bin/bash

. ../MasterTest.sh

CleanFiles avg.in output.dat avg.dat

INPUT="-i avg.in"
cat > avg.in <<EOF
readdata perres.peptide.dat
avg perres.peptide.dat out output.dat name V
avg perres.peptide.dat oversets out avg.dat name Over
runanalysis
EOF
RunCpptraj "Data set averaging and averaging over all sets test."
DoTest output.dat.save output.dat
DoTest avg.dat.save avg.dat

EndTest
exit 0
