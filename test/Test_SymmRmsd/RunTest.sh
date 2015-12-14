#!/bin/bash

. ../MasterTest.sh

CleanFiles rms.in rmsd?.dat rmsd.dat rotate.in

Rotate() {
  INPUT="-i rotate.in"
  cat > rotate.in <<EOF
parm AFV.parm7
trajin AFV.rst7
makestructure chi2:2:CA:CB:CG:CD1:0
#trajout AFV.rotate2.rst7
makestructure chi3:3:N:CA:CB:CG1:-120
#trajout AFV.rotate3.rst7
trajout AFV.rotate23.rst7
EOF
  RunCpptraj "Rotate res2 of AFV"
}

Rms() {
  INPUT="-i rms.in"
  cat > rms.in <<EOF
parm AFV.parm7
reference AFV.rst7
trajin AFV.rotate2.rst7
trajin AFV.rotate3.rst7 
trajin AFV.rotate23.rst7 
rmsd NoSymm !@H= reference out rmsd.dat
symmrmsd Symm !@H= reference out rmsd.dat
EOF
  RunCpptraj "SymmRmsd"
  DoTest rmsd.dat.save rmsd.dat
}

#Rotate
Rms

EndTest
exit 0
