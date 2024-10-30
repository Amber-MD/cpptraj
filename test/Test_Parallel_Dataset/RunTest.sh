#!/bin/bash

TESTNAME='Test handling of read-in vs generated data sets'

# The purpose of these tests is to make sure that CPPTRAJ consistently
# handles read-in data sets and generated data sets the same way in
# both serial and parallel.

. ../MasterTest.sh

CleanFiles cpptraj.in TZ2.RM.dat TZ2.rotate.dat TZ2.rotate.fromset.dat TZ2.rotate2.dat

INPUT='-i cpptraj.in'

UNITNAME="Handling of generated dataset with 'rotate' command"
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > cpptraj.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 100
reference ../tz2.nc lastframe [LAST]

# Fit on protein
rmsd TZ2 ref [LAST] :1-12@CA savematrices matricesout TZ2.RM.dat
# Extract rotations from matrices
rotate calcfrom TZ2[RM] name Rot out TZ2.rotate.dat
run

# Extract rotations from the generated set (set should count as read-in now)
rotate calcfrom TZ2[RM] name Rot2 out TZ2.rotate2.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest TZ2.RM.dat.save TZ2.RM.dat
  DoTest TZ2.rotate.dat.save TZ2.rotate.dat
  DoTest TZ2.rotate.dat.save TZ2.rotate2.dat -I \#Frame
fi

UNITNAME="Handling of read-in dataset with 'rotate' command"
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > cpptraj.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 100
readdata TZ2.RM.dat name TZ2RM mat3x3

# Fit on protein
#rmsd TZ2 first :1-12@CA
# Extract rotations from matrices
rotate calcfrom TZ2RM name Rot out TZ2.rotate.fromset.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest TZ2.rotate.dat.save TZ2.rotate.fromset.dat
fi

EndTest
