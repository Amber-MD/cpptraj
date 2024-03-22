#!/bin/bash

. ../MasterTest.sh

CleanFiles tica.in ticadebug.dat

INPUT='-i tica.in'

MASS=mass
cat > tica.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.crd name TZ2

set MASK = :1-3@N,CA,C

crdaction TZ2 matrix name M1 mwcovar out M1.dat \$MASK $MASS

runanalysis tica crdset TZ2 mask \$MASK lag 1 debugfile ticadebug.dat $MASS

list dataset
EOF
RunCpptraj "TICA test."

EndTest
exit 0
