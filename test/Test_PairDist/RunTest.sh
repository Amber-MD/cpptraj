#!/bin/bash

. ../MasterTest.sh

in=Pr.in
out=Pr.dat

CleanFiles $in $out ortho.dat truncoct.dat twomask.dat

INPUT="-i $in"

cat > $in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd

pairdist P0 out $out mask "*" delta 0.1 maxdist 36.25
EOF

RunCpptraj "PairDist Test."
DoTest ${out}.save $out

UNITNAME='Pairdist test, orthogonal imaging'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > $in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
pairdist P1 out ortho.dat mask :WAT@O delta 0.1 maxdist 20.0
EOF
  RunCpptraj "$UNITNAME."
  DoTest ortho.dat.save ortho.dat
fi

UNITNAME='Pairdist test, nonorthogonal imaging'
CheckFor netcdf maxthreads 2
if [ $? -eq 0 ] ; then
  cat > $in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
pairdist P2 out truncoct.dat mask :WAT@O delta 0.1 maxdist 20.0
EOF
  RunCpptraj "$UNITNAME."
  DoTest truncoct.dat.save truncoct.dat
fi

UNITNAME='Pairdist test, 2 masks'
cat > $in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd
pairdist P0 out twomask.dat mask @O mask2 @H delta 0.1 maxdist 36.25
EOF
RunCpptraj "$UNITNAME."
DoTest twomask.dat.save twomask.dat

EndTest

exit 0
