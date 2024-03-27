#!/bin/bash

. ../MasterTest.sh

CleanFiles tica.in ticadebug.dat M1.dat

INPUT='-i tica.in'

MASS=mass
cat > tica.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.crd name TZ2

set MASK = :1-3@N,CA,C

#crdaction TZ2 matrix name M1 mwcovar out M1.dat \$MASK $MASS

#runanalysis tica crdset TZ2 mask \$MASK lag 1 debugfile ticadebug.dat $MASS

crdaction TZ2 matrix name M1 mwcovar out M1.dat :1-3@N,CA :1-3@C,O $MASS

runanalysis tica crdset TZ2 mask :1-3@N,CA mask2 :1-3@C,O lag 1 debugfile ticadebug.dat $MASS

list dataset
EOF
RunCpptraj "TICA test."

EndTest
exit 0
