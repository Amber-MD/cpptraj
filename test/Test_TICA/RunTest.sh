#!/bin/bash

. ../MasterTest.sh

CleanFiles tica.in

INPUT='-i tica.in'

cat > tica.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.crd name TZ2

crdaction TZ2 matrix name M1 covar out M1.dat @CA

runanalysis tica crdset TZ2 mask @CA lag 1
EOF
RunCpptraj "TICA test."

EndTest
exit 0
