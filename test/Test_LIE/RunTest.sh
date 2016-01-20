#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles lie.in test.out LIE.out

# Check libraries
CheckNetcdf
CheckZlib
CheckBzlib
MaxThreads 2
if [[ $? -eq 0 ]] ; then
  cat > lie.in <<EOF
trajin test.mdcrd
lie LIE :RAL out LIE.out cutvdw 12 cutelec 12
EOF
  INPUT="lie.in"
  TOP="sol.top"
  RunCpptraj "LIE Test"
  DoTest LIE.out.save LIE.out
fi
EndTest

exit 0
