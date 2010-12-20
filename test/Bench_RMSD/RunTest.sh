#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in ptraj.in rmsd.dat time.dat Timing_Results.dat 


#for ((ATOMS=10; ATOMS < 4049; ATOMS += 200)) ; do
  # Run each test 3 times
#  for NTEST in test1 test2 test3 ; do
    #echo $ATOMS
    #TRAJ=../ptraj.run0.crd
    TRAJ=../run0.nc
    cat > ptraj.in <<EOF
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
trajin $TRAJ
rms first out rmsd.dat :1-268@CA 
rms first out rmsd.dat  
rms first out rmsd.dat :191-211@N,CA,C 
rms first out rmsd.dat :1-100 
EOF
    cp ptraj.in cpptraj.in
    echo "noprogress" >> cpptraj.in
    #CPPTRAJ=/home/droe/bin/cpptraj
    TOP=../ChainA-tip3p.parm7
    INPUT="cpptraj.in"
    TIME="time"
    ERROR="time.dat"
    echo $CPPTRAJ
    RunCpptraj "CPPTRAJ: Timing of RMSD command."
    
    CPPTRAJ=`which ptraj`
    echo $CPPTRAJ
    INPUT="ptraj.in"
    RunCpptraj "PTRAJ: Timing of RMSD command."

    #DoTest d1.dat.save d1.dat 
    #CheckTest
    # Discard first test results
    #if [[ $NTEST = "test1" ]] ; then
    #  CleanFiles $ERROR
    #fi
#  done
#done
mv Test_Results.dat Timing_Results.dat
#EndTest
cat $ERROR
exit 0
