#!/bin/bash

. ../MasterTest.sh

CleanFiles mcurve.in curves.agr Results.dat

INPUT="-i mcurve.in"
cat > mcurve.in <<EOF
readdata crv.1-3.life.agr name LC
#debug analysis 5
runanalysis multicurve set LC resultsout Results.dat \
  name Fit nexp 2 form mexp A0=1.0 A1=-1.0 tol 0.0001 maxit 50
list datasets
writedata curves.agr LC Fit*
EOF
RunCpptraj "Multicurve test"
DoTest curves.agr.save curves.agr
EndTest
exit 0
