#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles nastruct.in BP.*.dat BPstep.*.dat bases.pdb baseaxes.pdb basepairaxes.pdb \
           Helix.*.dat Param.pdb

# Test 2
MaxThreads 3 "NAstruct tests"
if [[ $? -eq 0 ]] ; then
  INPUT="-i nastruct.in"
  cat > nastruct.in <<EOF
parm ../adh026.3.pdb
trajin ../adh026.3.pdb 
nastruct naout adh026.dat
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
fi
EndTest

exit 0
