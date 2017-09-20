#!/bin/bash

. ../MasterTest.sh

CleanFiles create.in crd?.dat crd?.gnu crd?.rst7 crd.rst7

INPUT="-i create.in"
TESTNAME='COORDS creation tests'
Requires netcdf
# Test COORDS creation and processing
cat > create.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
atomicfluct out crd1.dat byatom bfactor
createcrd crd1
run
crdfluct crdset crd1 out crd1.dat window 20 bfactor
runanalysis
EOF
RunCpptraj "COORDS data set creation and CRDFLUCT test."
DoTest crd1.dat.save crd1.dat
 
# Test velocities
UNITNAME='COORDS data set creation with velocities test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > create.in <<EOF
parm ../rGACC.full.parm7
trajin ../rGACC.full.nc
trajout crd1.rst7 
trajout crd0.rst7 title " " nobox
createcrd crd1
run
crdout crd1 crd.rst7 title " " nobox time0 50961
EOF
  RunCpptraj "$UNITNAME"
  DoTest crd0.rst7 crd.rst7
  DoTest crd1.rst7.save crd1.rst7
fi

# Test TRAJ creation and processing
cat > create.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
loadtraj name crd1
atomicfluct out crd1.dat byatom bfactor
run
crdfluct crdset crd1 out crd1.dat window 20 bfactor
runanalysis
EOF
RunCpptraj "TRAJ data set creation and CRDFLUCT test."
DoTest crd1.dat.save crd1.dat

# Test TRAJ creation with velocities
cat > create.in <<EOF
parm ../rGACC.full.parm7
loadtraj ../rGACC.full.nc name crd2
crdout crd2 crd2.rst7 time0 50961
EOF
RunCpptraj "TRAJ data set creation with velocities test."
DoTest crd1.rst7.save crd2.rst7

EndTest

exit 0
