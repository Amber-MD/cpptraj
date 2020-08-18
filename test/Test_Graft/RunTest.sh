#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in Final.graft.mol2
TESTNAME='Graft test'
Requires notparallel

INPUT="-i cpptraj.in"
# Combine Tyr FF14SB backbone + CB with PRY fragment
cat > cpptraj.in <<EOF
set TYRFILE = ../Test_CombineCrd/Tyr.mol2
# Load the target
parm \$TYRFILE 
loadcrd \$TYRFILE parmindex 0 name TYR

# Load the source
set PRYFILE= ../Test_CombineCrd/PRY-gauss-fragment.mol2
parm \$PRYFILE
loadcrd \$PRYFILE parmindex 1 name PRY

graft \
  tgt TYR \
  src PRY \
  name Final \
  srcfitmask @O1,C5,C6,C4,H6,H5,C7,C3,H7,H4,C2 \
  tgtfitmask @OH,CZ,CE1,CE2,HE1,HE2,CD1,CD2,HD1,HD2,CG \
  srcmask !(@1-4) \
  tgtmask @C,O,CA,HA,N,H,CB,HB2,HB3 \
  bond :1@CB,:1@C2
crdout Final Final.graft.mol2

EOF
RunCpptraj "$TESTNAME"
DoTest Final.graft.mol2.save Final.graft.mol2

EndTest
exit 0
