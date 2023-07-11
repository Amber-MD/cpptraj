#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles nastruct.in BP.*.dat BPstep.*.dat bases.pdb baseaxes.pdb basepairaxes.pdb \
           Helix.*.dat Param.pdb SS.mol1.dat SS.mol1.selected.dat \
           axes.bases.pdb axes.bp.mol2 axes.step.dcd axes.step.parm7

# Test 2
TESTNAME='NAstruct tests'
Requires maxthreads 3

INPUT="-i nastruct.in"
cat > nastruct.in <<EOF
parm ../adh026.3.pdb
trajin ../adh026.3.pdb 
nastruct naout adh026.dat \
  axesout axes.bases.pdb \
  bpaxesout axes.bp.mol2 \
  stepaxesout axes.step.dcd stepaxesparmout axes.step.parm7
nastruct naout baseref.dat baseref Atomic_G.pdb.nastruct
nastruct naout groove.dat groovecalc 3dna
nastruct naout GuessBP.dat guessbp
EOF
RunCpptraj "NAstruct command test."
DoTest BP.adh026.dat.save BP.adh026.dat
DoTest BPstep.adh026.dat.save BPstep.adh026.dat
DoTest Helix.adh026.dat.save Helix.adh026.dat
DoTest BP.adh026.dat.save BP.baseref.dat
DoTest BPstep.adh026.dat.save BPstep.baseref.dat
DoTest Helix.adh026.dat.save Helix.baseref.dat
DoTest BPstep.groove.dat.save BPstep.groove.dat
DoTest BP.adh026.dat.save BP.GuessBP.dat
DoTest BPstep.adh026.dat.save BPstep.GuessBP.dat
DoTest Helix.adh026.dat.save Helix.GuessBP.dat
DoTest axes.bases.pdb.save axes.bases.pdb
DoTest axes.bp.mol2.save axes.bp.mol2
DoTest axes.step.parm7.save axes.step.parm7 -I %VERSION

# Single strand
cat > nastruct.in <<EOF
parm ../adh026.3.pdb
trajin ../adh026.3.pdb
strip ^2
nastruct naout mol1.dat sscalc
EOF
RunCpptraj "NAstruct command test, single strand"
DoTest SS.mol1.dat.save SS.mol1.dat

# Single strand, selected residues
cat > nastruct.in <<EOF
parm ../adh026.3.pdb
trajin ../adh026.3.pdb
strip ^2
nastruct naout mol1.selected.dat sscalc resrange 2,3,5,6,7
EOF
RunCpptraj "NAstruct command test, single strand, selected residues"
DoTest SS.mol1.selected.dat.save SS.mol1.selected.dat

EndTest

exit 0
