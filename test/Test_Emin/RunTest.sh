#!/bin/bash

. ../MasterTest.sh

TESTNAME='Emin tests'

CleanFiles emin.in cpptraj.ene.dat cpptraj.emin.nc \
           omm.ene.dat ommangle.ene.dat

INPUT='-i emin.in'

UNITNAME='Basic energy minimization test.'
cat > emin.in <<EOF
parm O2mol.parm7
loadcrd O2mol.rst7 name O2mol

emin crdset O2mol nsteps 100 out cpptraj.ene.dat trajoutname cpptraj.emin.nc
EOF
RunCpptraj "$UNITNAME"
DoTest cpptraj.ene.dat.save cpptraj.ene.dat

UNITNAME='OpenMM basic energy minimization tests'
CheckFor openmm
if [ $? -eq 0 ] ; then
  # Bonds
  cat > emin.in <<EOF
parm O2mol.parm7
loadcrd O2mol.rst7 name O2mol
emin crdset O2mol nsteps 100 out omm.ene.dat openmm
EOF
  RunCpptraj "$UNITNAME (bonds)"
  DoTest cpptraj.ene.dat.save omm.ene.dat
  # Angles
  cat > emin.in <<EOF
parm HOHmol.parm7
loadcrd HOHmol.rst7 name HOHmol
emin crdset HOHmol nsteps 100 out ommangle.ene.dat openmm
EOF
  RunCpptraj "$UNITNAME (angles)"
  DoTest ommangle.ene.dat.save ommangle.ene.dat
fi


EndTest
