#!/bin/bash

. ../MasterTest.sh

if [[ ! -e Restart ]] ; then
  mkdir Restart
fi

# Clean
CleanFiles lie.in test.out

# Check libraries
CheckNetcdf
CheckZlib
CheckBzlib

cat > lie.in <<EOF
trajin ../../mmpbsa_py/EstRAL_Files/test.mdcrd
lie LIE :RAL out LIE.out cutvdw 12 cutelec 12
EOF

INPUT="lie.in"
TOP="../../mmpbsa_py/EstRAL_Files/sol.top"
RunCpptraj "LIE Test"

DoTest LIE.out.save LIE.out
CheckTest

EndTest

exit 0
