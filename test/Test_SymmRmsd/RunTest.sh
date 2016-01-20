#!/bin/bash

. ../MasterTest.sh

CleanFiles rms.in rmsd?.dat rmsd.dat rotate.in srms.in *.rmsd.dat \
           TYR.remap.crd ASP.remap.crd GLU.remap.crd

MaxThreads 3 "Symmetry-corrected RMSD tests"
if [[ $? -ne 0 ]] ; then
  EndTest
  exit 0
fi

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

STest() {
  INPUT="-i srms.in"
  cat > srms.in <<EOF
parm $1.parm7
trajin $1.nc
rms first out $1.rmsd.dat NOSYMM
symmrmsd first out $1.rmsd.dat SYMM remap
trajout $1.remap.crd
EOF
  RunCpptraj "$1 symmetry-corrected rmsd test."
  DoTest $1.rmsd.dat.save $1.rmsd.dat
  DoTest $1.remap.crd.save $1.remap.crd
}

#Rotate
Rms
MaxThreads 2 "ASP and GLU symmetry tests"
if [[ $? -eq 0 ]] ; then
  STest ASP
  STest GLU
fi
STest TYR

EndTest
exit 0
