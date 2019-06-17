#!/bin/bash

. ../MasterTest.sh

CleanFiles dih.in dihrms.dat toref.dat previous.dat totraj.dat

TESTNAME='Dihedral RMSD tests'
Requires maxthreads 10

INPUT="-i dih.in"

UNITNAME='Dihedral RMSD tests (first, reference, reference traj.)'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
reference ../tz2.nc 1 [MyRef]
dihrms ToFirst out dihrms.dat noheader phi psi
dihrms ToRef ref [MyRef] out toref.dat noheader phi psi
# NOTE: All calculated values for same traj should be 0.0
dihrms ToTraj reftraj ../tz2.nc out totraj.dat phi psi
EOF
  RunCpptraj "$UNITNAME"
  DoTest dihrms.dat.save dihrms.dat
  DoTest dihrms.dat.save toref.dat
  DoTest totraj.dat.save totraj.dat
fi

UNITNAME='Dihedral RMSD test (to previous)'
CheckFor notparallel
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
dihrms ToPrev previous out previous.dat phi psi
EOF
  RunCpptraj "$UNITNAME"
  DoTest previous.dat.save previous.dat
fi

EndTest

exit 0
