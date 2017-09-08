#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in Final.PRY.mol2 combinedCrd.crd combinedParm.parm7 \
           E?.dat FabI.NDP.TCS.parm7
TESTNAME='Combine COORDS tests'
Requires notparallel

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

UNITNAME='Combine COORDS with box info, split/combine COORDS tests'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  # Combine triclosan from FtuFabI + NAD + TCS with trpzip2
  cat > cpptraj.in <<EOF
parm ../FtuFabI.NAD.TCL.parm7
loadcrd ../FtuFabI.NAD.TCL.nc name TCS
crdaction TCS strip !:TCS

parm ../tz2.truncoct.parm7
loadcrd ../tz2.truncoct.nc name TZ2 parmindex 1
crdaction TZ2 strip :WAT

combinecrd TZ2 TCS parmname combinedParm crdname combinedCrd
parmwrite out combinedParm.parm7 parm combinedParm
crdout combinedCrd combinedCrd.crd
EOF
  RunCpptraj "Combine COORDS test with box info."
  DoTest combinedParm.parm7.save combinedParm.parm7 -I %VERSION
  DoTest combinedCrd.crd.save combinedCrd.crd

  # Split FtuFabI+NDP+TCS apart, then put back together.
  CRD='../FtuFabI.NAD.TCL.nc 1 3'
  cat > cpptraj.in <<EOF
parm ../FtuFabI.NAD.TCL.parm7

loadcrd $CRD name FabI
# Get initial system energies
crdaction FabI energy out E1.dat noheader
crdaction FabI strip :NDP,TCS

loadcrd $CRD name NDP
crdaction NDP strip !:NDP

loadcrd $CRD name TCS
crdaction TCS strip !:TCS

combinecrd FabI NDP TCS parmname combinedParm crdname combinedCrd
# Get recombined system energies. Should match exactly.
crdaction combinedCrd energy out E2.dat noheader
parmwrite out FabI.NDP.TCS.parm7
EOF
  RunCpptraj "Split coords and recombine test."
  DoTest E1.dat E2.dat
  DoTest ../FtuFabI.NAD.TCL.parm7 FabI.NDP.TCS.parm7 -I %VERSION
fi

EndTest
exit 0
