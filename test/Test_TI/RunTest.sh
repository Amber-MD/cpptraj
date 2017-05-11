#!/bin/bash

. ../MasterTest.sh

CleanFiles ti.in skip.agr curve.agr bs.dat incr.agr incr.crv.agr bs.crv.agr \
           avg.dat.save avg.dat

INPUT='-i ti.in'

cat > ti.in <<EOF
readdata dvdl.dat index 1 name TI
ti TI nq 9 name Curve out skip.agr curveout curve.agr bsout bs.dat \
  nskip 0,5,10,15,20,30,40,50
ti TI nq 9 name Avg
ti TI nq 9 name Increment avgincrement 10 out incr.agr curveout incr.crv.agr
ti TI nq 9 name Bootstrap bs_samples 20 bs_seed 10 bs_pts 70 out bs.dat curveout bs.crv.agr
runanalysis
list dataset
printdata Avg
writedata avg.dat.save Curve[TIcurve]:0 noheader
writedata avg.dat      Avg[TIcurve]     noheader
EOF
RunCpptraj "TI analysis test."
DoTest skip.agr.save skip.agr
DoTest curve.agr.save curve.agr
DoTest avg.dat.save avg.dat
DoTest bs.dat.save bs.dat

EndTest
exit 0
