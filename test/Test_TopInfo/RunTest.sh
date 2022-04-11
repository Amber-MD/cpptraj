#!/bin/bash

. ../MasterTest.sh

CleanFiles info.in atoms.dat residues.dat bonds.dat angles.dat dihedrals.dat \
           molecules.dat masscharge.dat values.dat molshort.dat molselect.dat \
           molselect2.dat ChargeMass.dat

INPUT="-i info.in"
cat > info.in <<EOF
parm ../tz2.parm7
atoms :3
atoms :3 out atoms.dat

resinfo *
resinfo out residues.dat *
resinfo short *
resinfo short out residues.dat *

bonds :1
bonds :1 out bonds.dat
bonds @%N3 @%H
bonds @%N3 @%H out bonds.dat

angles @1
angles @1 out angles.dat
angles @%H @%N3 @%CT
angles @%H @%N3 @%CT out angles.dat

dihedralinfo @1
dihedralinfo @1 out dihedrals.dat
dihedralinfo @N @CA @CB @%H1
dihedralinfo @N @CA @CB @%H1 out dihedrals.dat

mass out masscharge.dat name Mass *
charge out masscharge.dat name Charge *
writedata ChargeMass.dat Charge Mass noheader noxcol

parm ../dna30.parm7
molinfo !:WAT 1
molinfo !:WAT out molecules.dat 1
molinfo short 1 *
molinfo short out molshort.dat 1 *
resinfo ^2,5-7,100 out molselect.dat parm dna30.parm7
atoms ^1:DC@P,O?P out molselect2.dat parmindex 1
quit
EOF
RunCpptraj "Topology info print test."
DoTest atoms.dat.save atoms.dat
DoTest residues.dat.save residues.dat
DoTest bonds.dat.save bonds.dat
DoTest angles.dat.save angles.dat
DoTest dihedrals.dat.save dihedrals.dat
DoTest masscharge.dat.save masscharge.dat
DoTest molecules.dat.save molecules.dat
DoTest molshort.dat.save molshort.dat
DoTest molselect.dat.save molselect.dat
DoTest molselect2.dat.save molselect2.dat
DoTest ChargeMass.dat.save ChargeMass.dat

UNITNAME='Topology info with reference coords test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > info.in <<EOF
parm ../tz2.parm7
reference ../tz2.nc 50
bonds @10 out values.dat reference
angles @10 out values.dat reference
dihedrals @10 out values.dat reference
EOF
  RunCpptraj "$UNITNAME"
  DoTest values.dat.save values.dat
fi

EndTest
exit 0
