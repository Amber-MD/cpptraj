#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles closest.in Closest.pdb closest2.in first.Closest.pdb \
           closestmols.dat closest.tz2.truncoct.parm7 imaged.pdb \
           first.Closest.rst7 all.Closest.pdb

# Test 1 - Closest, first slovent atom only
CheckNetcdf
cat > closest.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
closest 10 :1-13 first closestout closestmols.dat name CL outprefix closest
trajout first.Closest.pdb pdb nobox 
EOF
INPUT="-i closest.in"
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
INPUT="-i closest.in"
RunCpptraj "Closest command test using all solvent atoms."
DoTest all.Closest.pdb.save all.Closest.pdb 

## Check, wrap the closest waters.
#cat > closest.in <<EOF
#parm closest.tz2.truncoct.parm7
#trajin first.Closest.rst7
#center :1-13
#image
#trajout imaged.pdb pdb
#EOF
#RunCpptraj "Re-imaging"

# Long test, disabled for now.
#cat > closest2.in <<EOF
#noprogress
#parm ../tz2.truncoct.parm7
#trajin ../tz2.truncoct.nc 1 1
#closest 10 :1-13
#trajout Closest.pdb pdb
#EOF
#INPUT="-i closest2.in"
#RunCpptraj "Closest command test using all solvent atoms."
#DoTest Closest.pdb.save Closest.pdb

CheckTest

EndTest

exit 0
