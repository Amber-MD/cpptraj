#!/bin/bash

. ../MasterTest.sh

CleanFiles keep.in solvavg.dat uvseries.dat hb.dat b.dat \
           keep.parm7 keep.dcd keep.crd \
           keep.10.11.parm7 keep.10.11.crd \
           res1.tz2.crd \
           hb.2.dat keep.onlyone.crd \
           keep.two.crd

INPUT='keep.in'

TESTNAME='Keep command tests.'

UNITNAME='Keep bridging water test'
CheckFor netcdf maxthreads 1
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
CheckFor netcdf maxthreads 1
if [ $? -eq 0 ] ; then
  cat > keep.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc

readdata hb.dat

keep keepmask :10,11 bridgedata hb.dat:5 parmout keep.10.11.parm7 nobridgewarn
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

UNITNAME='Keep only 1 bridging water test'
CheckFor netcdf maxthreads 1
if [ $? -eq 0 ] ; then
  cat > keep.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
# First pass, generate bridge time series
hbond hb solventacceptor :WAT@O solventdonor :WAT out hb.2.dat
run
#writedata b.dat hb[bridge_10_11]
# Second pass, retain only frames where the bridge is present

keep bridgedata hb[ID] nbridge 1 
trajout keep.onlyone.crd
run
EOF
  RunCpptraj "$UNITNAME"
  DoTest keep.onlyone.crd.save keep.onlyone.crd
fi

UNITNAME='Keep 2 bridging waters at specified residues test'
CheckFor netcdf maxthreads 1
if [ $? -eq 0 ] ; then
  cat > keep.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
# First pass, generate bridge time series
hbond hb solventacceptor :WAT@O solventdonor :WAT
run
#writedata b.dat hb[bridge_10_11]
# Second pass, retain only frames where the bridge is present

keep bridgedata hb[ID] nbridge 2 bridgeresonly 10,11,3,5
trajout keep.two.crd
run
EOF
  RunCpptraj "$UNITNAME"
  DoTest keep.two.crd.save keep.two.crd
fi

EndTest
