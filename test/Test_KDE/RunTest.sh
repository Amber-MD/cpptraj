#!/bin/bash

. ../MasterTest.sh

CleanFiles kde.in kde.dat

INPUT="-i kde.in"
cat > kde.in <<EOF
readdata ../Test_SPAM/spampure.dat.save name SPAM
runanalysis kde SPAM min -37 max 5 bins 100 out kde.dat name KDE_SPAM
EOF
RunCpptraj "KDE test."
DoTest kde.dat.save kde.dat

EndTest
exit 0
