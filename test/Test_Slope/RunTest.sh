#!/bin/bash

. ../MasterTest.sh

CleanFiles curve.in curve.agr

INPUT="-i curve.in"

cat > curve.in <<EOF
readdata ../Test_Corr/corr.dat.save name Corr
slope Corr out curve.agr prec 16.8 name Fwd type forward
slope Corr out curve.agr prec 16.8 name Bck type backward
slope Corr out curve.agr prec 16.8 name Cnt type central
EOF
RunCpptraj "Finite difference test."
DoTest curve.agr.save curve.agr

EndTest
exit 0
