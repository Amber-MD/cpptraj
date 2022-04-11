#!/bin/bash

. ../MasterTest.sh

CleanFiles phipsi.in dihedrals.dat phipsi.dat

TESTNAME='Phipsi test'
Requires netcdf

if [ $? -eq 0 ] ; then
  INPUT="-i phipsi.in"
  cat > phipsi.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
multidihedral DPDP phi psi
run
phipsi DPDP[phi] DPDP[psi] out phipsi.dat resrange 1-22
writedata dihedrals.dat DPDP[phi]:2 DPDP[psi]:2
EOF
  RunCpptraj "$TESTNAME"
  DoTest phipsi.dat.save phipsi.dat
fi

EndTest
exit 0
