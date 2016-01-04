#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles closest.in first.Closest.pdb closestmols.dat \
           closest.tz2.truncoct.parm7 all.Closest.pdb closest10.center2_4.crd

INPUT="-i closest.in"
CheckNetcdf
if [[ -z $DO_PARALLEL ]] ; then
  # Test 1 - Closest, first solvent atom only
  cat > closest.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
closest 10 :1-13 first closestout closestmols.dat name CL outprefix closest
trajout first.Closest.pdb pdb nobox 
EOF
  RunCpptraj "Closest command test using first solvent atom."
  DoTest first.Closest.pdb.save first.Closest.pdb
  DoTest closestmols.dat.save closestmols.dat
  # Tell diff to ignore the VERSION line
  DoTest closest.tz2.truncoct.parm7.save closest.tz2.truncoct.parm7 -I %VERSION

  # Test 2 - Closest, all solvent atoms
  cat > closest.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
closest 10 :1-13
trajout all.Closest.pdb pdb nobox
EOF
  RunCpptraj "Closest command test using all solvent atoms."
  DoTest all.Closest.pdb.save all.Closest.pdb 
fi

# Test 3 - Closest atoms to mask center
cat > closest.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
closest 10 :2,4 center
trajout closest10.center2_4.crd nobox
EOF
RunCpptraj "Closest command test, using mask center"
DoTest closest10.center2_4.crd.save closest10.center2_4.crd

EndTest

exit 0
