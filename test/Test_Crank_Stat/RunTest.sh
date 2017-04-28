#!/bin/bash

. ../MasterTest.sh
CleanFiles crank.in crank.dat results.dat stat.dat
CheckNetcdf
INPUT="crank.in"
TOP=../tz2.parm7
cat > crank.in <<EOF
trajin ../tz2.nc
dihedral phi2 :1@C :2@N :2@CA :2@C type phi
dihedral phi3 :2@C :3@N :3@CA :3@C type phi
analyze crankshaft angle phi2 phi3 info "residue 2-3 phi/phi" out crank.dat results results.dat
analyze statistics all out stat.dat
EOF
RunCpptraj "Crankshaft test"
DoTest crank.dat.save crank.dat
DoTest results.dat.save results.dat
DoTest stat.dat.save stat.dat

EndTest
exit 0
