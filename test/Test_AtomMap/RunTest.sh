#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles atommap.in initial.mol2 atommap.dat reordered.pdb reordered.mol2 \
           fit.mol2 rmsd.dat map.chm_to_amb.dat mapped.pdb.? rmsout.dat \
           map.byres.chm_to_amb.dat

RequiresMaxThreads 3 "Atom map tests"

INPUT="-i atommap.in"
MaxThreads 1 "AtomMap Test"
if [ $? -ne 0 ] ; then
  SkipCheck "AtomMap Test"
else
  # Test 1
  cat > atommap.in <<EOF
noprogress
parm xtallig.mol2
reference xtallig.mol2

parm start.mol2
reference start.mol2 parmindex 1

outtraj initial.mol2 mol2

# xtallig -> start
atommap xtallig.mol2 start.mol2 mapout atommap.dat

trajin xtallig.mol2
outtraj reordered.pdb pdb chainid X model
outtraj reordered.mol2 mol2
rms 1g9v refindex 1 out rmsd.dat
trajout fit.mol2 mol2
EOF
  RunCpptraj "AtomMap Test"
  DoTest initial.mol2.save initial.mol2
  DoTest atommap.dat.save atommap.dat
  DoTest reordered.pdb.save reordered.pdb
  DoTest reordered.mol2.save reordered.mol2
  DoTest rmsd.dat.save rmsd.dat
  DoTest fit.mol2.save fit.mol2
fi

# Test 2
MaxThreads 2 "AtomMap with 'rmsfit'"
if [ $? -ne 0 ] ; then
  SkipCheck "AtomMap Test with 'rmsfit'"
else
  cat > atommap.in <<EOF
parm xtallig.mol2
reference xtallig.mol2

parm start.mol2
reference start.mol2 parm start.mol2

atommap xtallig.mol2 start.mol2 rmsfit rmsout rmsout.dat 1g9v

trajin xtallig.mol2
EOF
  RunCpptraj "AtomMap with 'rmsfit'"
  DoTest rmsd.dat.save rmsout.dat
fi

# Test 3
MaxThreads 3 "Atom map charmm->amber atom order"
if [ $? -ne 0 ] ; then
  SkipCheck "Atom map charmm->amber atom order"
else
  cat > atommap.in <<EOF
parm cg-amb.topo
reference cg-amb.crds
parm cg-chm.topo
reference cg-chm.crds parmindex 1
atommap cg-amb.crds cg-chm.crds mapout map.chm_to_amb.dat
trajin cg.crd
trajout mapped.pdb pdb multi
EOF
  RunCpptraj "Atom map charmm->amber atom order"
  DoTest map.chm_to_amb.dat.save map.chm_to_amb.dat

  # Test 4
  cat > atommap.in <<EOF
parm cg-amb.topo
reference cg-amb.crds
parm cg-chm.topo
reference cg-chm.crds parmindex 1
atommap cg-amb.crds cg-chm.crds mapout map.byres.chm_to_amb.dat maponly mode byres
EOF
  RunCpptraj "Atom map charmm->amber atom order, by residue"
  DoTest map.chm_to_amb.dat.save map.byres.chm_to_amb.dat
fi
EndTest

exit 0
