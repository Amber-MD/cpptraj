#!/bin/bash

. ../MasterTest.sh

CleanFiles ti.in skip.agr curve.agr bs.dat

INPUT='-i ti.in'

cat > ti.in <<EOF
readdata dvdl.dat index 1 name TI
ti TI nq 9 name Curve out skip.agr curveout curve.agr bsout bs.dat \
  nskip 0,5,10,15,20,30,40,50 \
  bs_samples 20
runanalysis
EOF
RunCpptraj "TI analysis test."
DoTest skip.agr.save skip.agr
DoTest curve.agr.save curve.agr
DoTest bs.dat.save bs.dat

EndTest
exit 0
