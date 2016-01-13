#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles jcoupling.in Jcoupling.dat jc.dat

CheckNetcdf
INPUT="-i jcoupling.in"
# Test 1
MaxThreads 1 "J-Coupling single frame test."
if [[ $? -eq 0 ]] ; then
  cat > jcoupling.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc 1 1
jcoupling outfile Jcoupling.dat kfile Karplus.txt
EOF
  RunCpptraj "J-Coupling command test."
  DoTest Jcoupling.dat.save Jcoupling.dat
fi

cat > jcoupling.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
jcoupling out jc.dat kfile Karplus.txt :2
EOF
RunCpptraj "J-Coupling extended test."
DoTest jc.dat.save jc.dat

EndTest

exit 0
