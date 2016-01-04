#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles general.in *.rst7.? framerange.in framerange.crd
NotParallel "Trajout Frame Range"
if [[ $? -eq 1 ]] ; then
  echo ""
  exit 0
fi
CheckNetcdf
# Test 1
cat > general.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc 
trajout test.rst7 restart onlyframes 3,5,6-8 time0 1.0
EOF
INPUT="-i general.in"
RunCpptraj "Trajout Frame Range [Amber Restart]"
DoTest test.rst7.3.save test.rst7.3 
DoTest test.rst7.5.save test.rst7.5 
DoTest test.rst7.6.save test.rst7.6 
DoTest test.rst7.7.save test.rst7.7 
DoTest test.rst7.8.save test.rst7.8 

# Test 2
cat > framerange.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc 
trajout framerange.crd onlyframes 3,5,6-8
EOF
INPUT="-i framerange.in"
RunCpptraj "Trajout Frame Range [Amber Coordinates]"
DoTest framerange.crd.save framerange.crd

CheckTest
EndTest

exit 0
