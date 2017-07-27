#!/bin/bash

. ../MasterTest.sh

CleanFiles order.in ol?.dat ol3.agr dopc.agr

INPUT="-i order.in"

MaxThreads 1 "New Lipid Order Parameter Test."
if [ $? -eq 0 ] ; then
  cat > order.in <<EOF
# crd/top courtesy of Callum Dickson, Imperial College London
parm ../DOPC.parm7
trajin ../DOPC.rst7

lipidscd :OL,PC  out ol1.dat
#debug actions 3
lipidscd :OL2,PC out ol2.dat
lipidscd out ol3.agr xydy
precision ol1.dat 10 7
precision ol2.dat 10 7
precision ol3.agr 10 7
datafile ol1.dat xprec 3.0
datafile ol2.dat xprec 3.0
datafile ol3.agr xprec 3.0
EOF
  RunCpptraj "New Lipid Order Parameter Test."
  DoTest ol1.dat.save ol1.dat
  DoTest ol2.dat.save ol2.dat
  DoTest ol3.agr.save ol3.agr
fi

cat > order.in <<EOF
parm dopc.charmm.mol2
trajin dopc.charmm.nc
lipidscd DOPC out dopc.agr xydy prec 10.7 xprec 3.0
EOF
RunCpptraj "New Lipid Order Parameter Test, DOPC."

EndTest

exit 0
