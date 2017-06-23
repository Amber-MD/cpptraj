#!/bin/bash

. ../MasterTest.sh

CleanFiles log.in trepidx.agr

INPUT="-i log.in"
cat > log.in <<EOF
readdata trem.log crdidx 1,2,3,4,5,6,7,8
runanalysis remlog trem.log out trepidx.agr repidx name Tcrd

EOF
RunCpptraj "Replica log read/analyze test."
DoTest trepidx.agr.save trepidx.agr

EndTest
exit 0
