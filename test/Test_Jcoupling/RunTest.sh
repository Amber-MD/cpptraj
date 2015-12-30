#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles jcoupling.in Jcoupling.dat jc.dat

INPUT="-i jcoupling.in"
# Test 1
CheckNetcdf
cat > jcoupling.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc 1 1
jcoupling outfile Jcoupling.dat kfile Karplus.txt
EOF
RunCpptraj "J-Coupling command test."
DoTest Jcoupling.dat.save Jcoupling.dat

cat > jcoupling.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
jcoupling out jc.dat kfile Karplus.txt :2
EOF
RunCpptraj "J-Coupling extended test."
DoTest jc.dat.save jc.dat

EndTest

exit 0
