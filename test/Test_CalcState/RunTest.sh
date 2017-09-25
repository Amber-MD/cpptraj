#!/bin/bash

. ../MasterTest.sh

CleanFiles calcstate.in state.dat curve.dat stateout states.dat trans.dat

INPUT='-i calcstate.in'

cat > calcstate.in <<EOF
readdata Phi9.dat name DIH
calcstate name State \
  state C0,DIH:2,-105.0,-75.0 \
  state C1,DIH:2,-135.0,-120.0 \
  out state.dat \
  curveout curve.dat \
  stateout states.dat \
  transout trans.dat
EOF
RunCpptraj "Calcstate test"
DoTest state.dat.save state.dat
DoTest curve.dat.save curve.dat
DoTest states.dat.save states.dat
DoTest trans.dat.save trans.dat

EndTest
exit 0
