#!/bin/bash

. ../MasterTest.sh

CleanFiles keep.in solvavg.dat uvseries.dat hb.dat b.dat \
           keep.parm7 keep.dcd keep.crd \
           keep.10.11.parm7 keep.10.11.crd \
           res1.tz2.crd

INPUT='keep.in'

TESTNAME='Keep command tests.'

UNITNAME='Keep bridging water test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > keep.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
# First pass, generate bridge time series
hbond hb :10,11 solventacceptor :WAT@O solventdonor :WAT \
      solvout solvavg.dat bridgeout solvavg.dat bseries \
      series uvseries uvseries.dat out hb.dat
run
#writedata b.dat hb[bridge_10_11]
# Second pass, retain only frames where the bridge is present

keep bridgedata hb[ID] parmout keep.parm7
trajout keep.crd
run
EOF
  RunCpptraj "$UNITNAME"
  DoTest keep.crd.save keep.crd
fi

UNITNAME='Keep bridging water and selected residues test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > keep.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
hbond hb :10,11 solventacceptor :WAT@O solventdonor :WAT
run
keep keepmask :10,11 bridgedata hb[ID] parmout keep.10.11.parm7
trajout keep.10.11.crd
run
EOF
  RunCpptraj "$UNITNAME"
fi

UNITNAME='Basic keep test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > keep.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc
keep keepmask :1 nobox
trajout res1.tz2.crd
EOF
  RunCpptraj "$UNITNAME"
  DoTest ../Test_Strip/res1.tz2.crd.save res1.tz2.crd
fi

EndTest
