#!/bin/bash

. ../MasterTest.sh

CleanFiles log.in trepidx.agr mremdreptime.dat ph.repidx.agr ph.stats.dat

INPUT="-i log.in"
cat > log.in <<EOF
readdata trem.log crdidx 1,2,3,4,5,6,7,8
runanalysis remlog trem.log out trepidx.agr repidx name Tcrd

readdata rem.log.1.save rem.log.2.save dimfile remd.dim as remlog nosearch
remlog rem.log.1.save stats reptime mremdreptime.dat

readdata ph.rem.log name PH
runanalysis remlog PH out ph.repidx.agr repidx stats statsout ph.stats.dat \
                   reptime ph.stats.dat name pHrem


EOF
RunCpptraj "Replica log read/analyze test."
DoTest trepidx.agr.save trepidx.agr
DoTest mremdreptime.dat.save mremdreptime.dat
DoTest ph.repidx.agr.save ph.repidx.agr
DoTest ph.stats.dat.save ph.stats.dat


EndTest
exit 0
