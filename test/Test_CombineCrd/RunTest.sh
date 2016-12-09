#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in Final.PRY.mol2

INPUT="-i cpptraj.in"
# Combine Tyr FF14SB backbone + CB with PRY fragment
cat > cpptraj.in <<EOF
parm Tyr.mol2
reference Tyr.mol2 parm Tyr.mol2

# Load the fragment, fit it on top of the existing tyrosine atoms
parm PRY-gauss-fragment.mol2 [pryparm]
loadcrd PRY-gauss-fragment.mol2 parm [pryparm] PRY
crdaction PRY strip @1-4
#crdaction PRY center @C2 origin
crdaction PRY rms reference @O1,C5,C6,C4,H6,H5,C7,C3,H7,H4,C2 @OH,CZ,CE1,CE2,HE1,HE2,CD1,CD2,HD1,HD2,CG

loadcrd Tyr.mol2 parm Tyr.mol2 TYR
crdaction TYR strip !@C,O,CA,HA,N,H,CB,HB2,HB3
#crdaction TYR center @CB origin
#crdaction TYR translate x 1.0 z -1.0

combinecrd TYR PRY parmname Final.PRY crdname Final
crdout Final Final.PRY.mol2
EOF

RunCpptraj "Combine COORDS test."
DoTest Final.PRY.mol2.save Final.PRY.mol2
EndTest
exit 0
