#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles closest.in Closest.pdb closest2.in first.Closest.pdb \
           closestmols.dat closest.tz2.truncoct.parm7 imaged.pdb \
           first.Closest.rst7 all.Closest.pdb center.closest.pdb

# Test 1 - Closest, first slovent atom only
CheckNetcdf
INPUT="-i closest.in"
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

# Test 3 - Closest atoms to mask center
cat > closest.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 5 5
closest 10 :2,4 center
trajout center.closest.pdb
EOF
RunCpptraj "Closest command test using mask center."
DoTest center.closest.pdb.save center.closest.pdb

CheckTest

EndTest

exit 0
