#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles solvent.in closest10.mol2

INPUT="-i solvent.in"
TESTNAME='Solvent test'
Requires maxthreads 1

# Mark methanol as solvent and use closest command 
cat > solvent.in <<EOF
parm AlaDipeptide.MEOH.parm7
solvent :MOH
trajin AlaDipeptide.MEOH.rst7
closest 10 !:MOH noimage
trajout closest10.mol2
EOF
RunCpptraj "Test solvent command with closest"
DoTest closest10.mol2.save closest10.mol2

EndTest

exit 0
