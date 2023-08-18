#!/bin/bash

. ../MasterTest.sh

CleanFiles tica.in

INPUT='-i tica.in'

cat > tica.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.crd name TZ2

set MASK = :1-3@CA

crdaction TZ2 matrix name M1 mwcovar out M1.dat \$MASK

runanalysis tica crdset TZ2 mask \$MASK lag 1 mass

list dataset
EOF
RunCpptraj "TICA test."

EndTest
exit 0
