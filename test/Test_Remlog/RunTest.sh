#!/bin/bash

. ../MasterTest.sh

CleanFiles log.in trepidx.agr mremdreptime.dat

INPUT="-i log.in"
cat > log.in <<EOF
readdata trem.log crdidx 1,2,3,4,5,6,7,8
runanalysis remlog trem.log out trepidx.agr repidx name Tcrd

readdata rem.log.1.save rem.log.2.save dimfile remd.dim as remlog nosearch
remlog rem.log.1.save stats reptime mremdreptime.dat


EOF
RunCpptraj "Replica log read/analyze test."
DoTest trepidx.agr.save trepidx.agr
DoTest mremdreptime.dat.save mremdreptime.dat

EndTest
exit 0
