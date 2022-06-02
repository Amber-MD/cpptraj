#!/bin/bash

# Using the 'prepareforleap' command to quickly parameterize glycoforms
# of the NIST monoclonal antibody.

rm Prepared.pdb *.in *.mol2 *.out *.parm7 *.rst7 *.log
if [ "$1" = 'clean' ] ; then
  exit 0
fi

CPPTRAJ=`which cpptraj`
if [ -z "$CPPTRAJ" ] ; then
  echo "CPPTRAJ not found."
  exit 1
fi

# See if LEaP is available.
LEAP=`which tleap`
if [ ! -z "$LEAP" ] ; then
  cat > leap.ff.in <<EOF
source leaprc.protein.ff14SB
source leaprc.GLYCAM_06j-1
EOF
  RUNLEAP='runleap leap.ff.in'
else
  echo "LEaP not present, will not run LEaP"
  RUNLEAP=''
fi

$CPPTRAJ <<EOF
parm mAb8671_2021_0526.pdb
loadcrd mAb8671_2021_0526.pdb name MyCrd

# Mutate PCA to GLN
change crdset MyCrd atomname from :PCA@OE to OE1
change crdset MyCrd resname from :PCA to GLN

prepareforleap \
  crdset MyCrd \
  name Prepared \
  pdbout Prepared.pdb \
  leapunitname Mol \
  out Prepared.leap.in $RUNLEAP \
  nowat noh

crdout Prepared Prepared.mol2 
EOF

