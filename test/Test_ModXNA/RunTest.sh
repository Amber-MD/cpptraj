#!/bin/bash

. ../MasterTest.sh

TESTNAME='ModXNA tests'

# These tests will help ensure compatability with the modxna.sh driver script.

CleanFiles cpptraj.in tmp.bonds tmp.Nucleotide.mol2

INPUT='-i cpptraj.in'

UNITNAME="Ensure 'sequence' works for ModXNA"
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > cpptraj.in <<EOF
parm tmp.BackboneSugar.mol2
# Check that the RK/REQ column that ModXNA expects is written
bonds @C1' @/H out tmp.bonds
loadcrd tmp.BackboneSugar.mol2 name BackboneSugar parm tmp.BackboneSugar.mol2
parm tmp.base-striped.mol2
loadcrd tmp.base-striped.mol2 name Base parm tmp.base-striped.mol2
dataset connect BackboneSugar headmask @C1'
dataset connect Base tailmask @N9
sequence Base BackboneSugar name Nucleotide
change crdset Nucleotide mergeres firstres 1 lastres 2
change crdset Nucleotide resname from * to AD3
change crdset Nucleotide oresnums of :1 min 1 max 1
# Charge on all atoms but O3' and HO3'
set Q1 = crdset Nucleotide charge inmask !@O3',HO3'
# Charge on all atoms plus O3'
Q2 = \$Q1 + -0.4
# What does charge on HO3' need to be to get target total 3cap charge?
Q3 = -0.679652 - \$Q2
change crdset Nucleotide charge of @HO3' to \$Q3
change crdset Nucleotide charge of @O3' to -0.4
crdout Nucleotide tmp.Nucleotide.mol2
EOF
  RunCpptraj "$UNITNAME"
  DoTest tmp.Nucleotide.mol2.save tmp.Nucleotide.mol2
  DoTest tmp.bonds.save tmp.bonds
fi

EndTest
