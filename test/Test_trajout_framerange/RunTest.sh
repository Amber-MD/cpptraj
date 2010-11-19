#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles *.rst7.? 
#rm distance.dat rmsd.dat rmsda.dat phi2.dat PhiPsi.dat test.crd a1.dat Restart/* test.nc
#rm r4.dat a2.dat.gz r2.dat r3-nofit.dat

cat > general.in <<EOF
noprogress
parm ../trpcage.parm7
trajin ../trpcage.nc 
trajout test.rst7 restart 3,5,6-8
EOF

INPUT="-i general.in"
#CPPTRAJ=/home/droe/bin/cpptraj
RunCpptraj "Trajout Frame Range"

DoTest test.rst7.2.save test.rst7.2 
DoTest test.rst7.4.save test.rst7.4 
DoTest test.rst7.5.save test.rst7.5 
DoTest test.rst7.6.save test.rst7.6 
DoTest test.rst7.7.save test.rst7.7 
CheckTest

EndTest

exit 0
