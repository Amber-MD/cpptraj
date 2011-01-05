#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles ptraj.in Closest.pdb Closest.pdb.1 cpptraj.in time.dat Timing_Results.dat 

# Test 1
CheckNetcdf
#cat > closest.in <<EOF
#noprogress
#parm ../ChainA-tip3p.parm7
#trajin ../run0.nc 1 1
#closest 10 :1-268 first
#trajout first.Closest.pdb pdb
#EOF
#INPUT="-i closest.in"
#RunCpptraj "Closest command test using first solvent atom."
#DoTest first.Closest.pdb.save first.Closest.pdb

cat > ptraj.in <<EOF
trajin ../run0.nc 1 1
closest 10 :1-268
trajout Closest.pdb pdb
EOF
cp ptraj.in cpptraj.in
echo "noprogress" >> cpptraj.in
TOP="../ChainA-tip3p.parm7"
INPUT="cpptraj.in"
TIME="time"
ERROR="time.dat"
echo $CPPTRAJ
RunCpptraj "Timing of Closest command using all solvent atoms."
#DoTest Closest.pdb.save Closest.pdb

CPPTRAJ=`which ptraj`
echo $CPPTRAJ
INPUT="ptraj.in"
RunCpptraj "PTRAJ: Timing of Closest command using all solvent atoms."

#CheckTest

#EndTest
mv Test_Results.dat Timing_Results.dat
cat $ERROR

exit 0
