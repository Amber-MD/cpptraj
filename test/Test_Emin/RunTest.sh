#!/bin/bash

. ../MasterTest.sh

TESTNAME='Emin tests'

CleanFiles emin.in cpptraj.ene.dat cpptraj.emin.nc \
           omm.ene.dat ommangle.ene.dat ommdihedral.ene.dat \
           ommnb.ene.dat omm.tz2.ene.dat omm.tz2.pme.dat \
           cpptraj.angle.dat cpptraj.dihedral.dat \
           cpptraj.nonbond.dat

INPUT='-i emin.in'

UNITNAME='Basic energy minimization test.'
cat > emin.in <<EOF
parm O2mol.parm7
loadcrd O2mol.rst7 name O2mol

emin crdset O2mol nsteps 100 out cpptraj.ene.dat #trajoutname cpptraj.emin.nc
EOF
RunCpptraj "$UNITNAME"
DoTest cpptraj.ene.dat.save cpptraj.ene.dat

UNITNAME='Energy minimization with angle term test.'
cat > emin.in <<EOF
parm HOHmol.parm7
loadcrd HOHmol.rst7 name HOHmol

emin crdset HOHmol nsteps 100 out cpptraj.angle.dat
EOF
RunCpptraj "$UNITNAME"
DoTest cpptraj.angle.dat.save cpptraj.angle.dat

UNITNAME='Energy minimization with dihedral and 1-4 terms test.'
cat > emin.in <<EOF
parm ETHmol.parm7
loadcrd ETHmol.rst7 name ETHmol

emin crdset ETHmol nsteps 100 out cpptraj.dihedral.dat
EOF
RunCpptraj "$UNITNAME"
DoTest cpptraj.dihedral.dat.save cpptraj.dihedral.dat

UNITNAME='Energy minimization with nonbond term test.'
# There are some small numerical differences on windows
# and with PGI compilers.
nbtestargs='-a 0.002'
cat > emin.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.rst7 1 1 name TZ2
emin crdset TZ2 nsteps 100 out cpptraj.nonbond.dat nonbond cut 9999
EOF
RunCpptraj "$UNITNAME"
DoTest cpptraj.nonbond.dat.save cpptraj.nonbond.dat $nbtestargs

# OpenMM tests
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
  # Dihedrals
  cat > emin.in <<EOF
parm ETHmol.parm7
loadcrd ETHmol.rst7 name ETHmol
emin crdset ETHmol nsteps 100 out ommdihedral.ene.dat openmm
EOF
  RunCpptraj "$UNITNAME (dihedrals)"
  DoTest ommdihedral.ene.dat.save ommdihedral.ene.dat
  # Nonbonds
  cat > emin.in <<EOF
parm O2box.parm7 
loadcrd O2box.rst7 1 1 name O2BOX
emin crdset O2BOX nsteps 100 out ommnb.ene.dat openmm
EOF
  RunCpptraj "$UNITNAME (nonbonds)"
  DoTest ommnb.ene.dat.save ommnb.ene.dat
  # All, non-periodic
  cat > emin.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.rst7 1 1 name TZ2
emin crdset TZ2 nsteps 100 out omm.tz2.ene.dat openmm
EOF
  RunCpptraj "$UNITNAME (Tz2)"
  DoTest omm.tz2.ene.dat.save omm.tz2.ene.dat
  # All, PME
  UNITNAME="$UNITNAME (Tz2 PME)"
  CheckFor netcdf
  if [ $? -eq 0 ] ; then
    cat > emin.in <<EOF
parm ../tz2.ortho.parm7
loadcrd ../tz2.ortho.nc 1 1 name TZ2
emin crdset TZ2 nsteps 100 out omm.tz2.pme.dat openmm !:WAT
EOF
    RunCpptraj "$UNITNAME"
    DoTest omm.tz2.pme.dat.save omm.tz2.pme.dat
  fi
fi

EndTest
