#!/bin/bash

. ../MasterTest.sh

CleanFiles cphstats.in

INPUT='-i cphstats.in'
TESTNAME='Constant pH stats / data set sort test'

UNITNAME='Ensemble data read / sort'
cat > cphstats.in <<EOF
readensembledata cpout.001 cpin cpin name PH
#readensembledata cpout.001 filenames cpout.002,cpout.003,cpout.004,cpout.005,cpout.006 name PH
list dataset
sortensembledata PH
list dataset
EOF
RunCpptraj "$UNITNAME"

EndTest
exit 0
