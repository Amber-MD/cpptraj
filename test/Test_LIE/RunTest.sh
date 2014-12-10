#!/bin/bash

. ../MasterTest.sh

if [[ ! -e Restart ]] ; then
  mkdir Restart
fi

# Clean
CleanFiles lie.in test.out LIE.out

# Check libraries
CheckNetcdf
CheckZlib
CheckBzlib

cat > lie.in <<EOF
trajin test.mdcrd
lie LIE :RAL out LIE.out cutvdw 12 cutelec 12
EOF

INPUT="lie.in"
TOP="sol.top"
RunCpptraj "LIE Test"

DoTest LIE.out.save LIE.out
CheckTest

EndTest

exit 0
