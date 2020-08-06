#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles jcoupling.in Jcoupling.dat jc.dat

TESTNAME='J-coupling tests'
Requires netcdf
INPUT="-i jcoupling.in"
# Test 1
UNITNAME='J-Coupling single frame test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > jcoupling.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc 1 1
jcoupling outfile Jcoupling.dat kfile ../../dat/Karplus.txt
EOF
  RunCpptraj "$UNITNAME"
  DoTest Jcoupling.dat.save Jcoupling.dat
fi

cat > jcoupling.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
jcoupling out jc.dat :2
EOF
RunCpptraj "J-Coupling extended test."
DoTest jc.dat.save jc.dat

EndTest

exit 0
