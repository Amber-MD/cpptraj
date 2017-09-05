#!/bin/bash

. ../MasterTest.sh

CleanFiles kde.in kde.dat kl.dat final.dat

INPUT="-i kde.in"
cat > kde.in <<EOF
readdata ../Test_SPAM/spampure.dat.save name SPAM index 1
runanalysis kde SPAM min -37 max 5 bins 100 out kde.dat name KDE_SPAM
EOF
RunCpptraj "KDE test."
DoTest kde.dat.save kde.dat

cat > kde.in <<EOF
parm ../DPDP.parm7
loadcrd ../DPDP.nc name DPDP
crdaction DPDP multidihedral PHI1 phi resrange 3 crdframes 1,50
crdaction DPDP multidihedral PHI2 phi resrange 3 crdframes 51,100
runanalysis kde PHI1[*] kldiv PHI2[*] klout kl.dat \
                min -133 max -52 step 5 out final.dat name KDE_Phi
EOF
RunCpptraj "KDE with Kullback-Leibler divergence test."
DoTest kl.dat.save kl.dat
DoTest final.dat.save final.dat

EndTest
exit 0
